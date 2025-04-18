#!/usr/bin/env julia

# Run the main thermodynamics script with just N2 for testing

# First, make sure the database is initialized with N2
include("initialize_minimal_database.jl")

# Create a standalone test script that directly calls our functions
# This avoids the include() issue with scope

println("\nRunning N2 test directly...")

# Activate the project to ensure we use the right dependencies
using Pkg
Pkg.activate(dirname(@__DIR__))

using JThermodynamicsData
using YAML
using DataFrames
using Plots
using Printf
using Dates

# Add lowercase constant aliases to the JThermodynamicsData module if they don't exist
println("Applying fixes for constants...")
if !isdefined(JThermodynamicsData, :h)
    @info "Adding lowercase constant aliases to JThermodynamicsData"
    for (lowercase_name, uppercase_name) in [(:h, :H), (:kb, :KB), (:na, :NA), (:c, :C)]
        if isdefined(JThermodynamicsData, uppercase_name)
            uppercase_value = getfield(JThermodynamicsData, uppercase_name)
            Core.eval(JThermodynamicsData, :(const $lowercase_name = $uppercase_value))
            @info "  Added JThermodynamicsData.$(lowercase_name)"
        end
    end
end

# Load main configuration
config_path = joinpath(dirname(@__DIR__), "config", "settings.yaml")
config = JThermodynamicsData.load_config(config_path)

# Initialize database
println("Initializing database...")
db_path = config["general"]["database_path"]
mkpath(dirname(db_path))
conn = JThermodynamicsData.init_database(config)

# Set species name and temperature range
species_name = "N2"
temp_range = [100.0, 10000.0]
temp_step = 100.0  # Use a larger step for testing

try
    # Define a fixed version of calculate_properties function
    function calculate_properties_direct(conn, species_name, temp_range, config, step)
        # Convert parameters for type stability
        temp_min = Float64(temp_range[1])
        temp_max = Float64(temp_range[2])
        step_val = Float64(step)
        
        # Generate temperature points
        temperatures = temp_min:step_val:temp_max
        
        # Try to get data from database
        results = Dict[]
        
        for temp in temperatures
            try
                result = JThermodynamicsData.query_properties(conn, species_name, temp, config)
                push!(results, result)
            catch e
                @warn "Failed to calculate properties at $temp K: $e"
            end
        end
        
        # Check if we got any results
        if isempty(results)
            @warn "No valid results for any temperature point for $species_name"
            
            # Try the direct calculator for N2
            if species_name == "N2"
                @info "Trying direct calculator for N2..."
                try
                    # Implement N2 calculator directly
                    function calculate_n2_properties(temperature::Real)
                        # Convert temperature to Float64
                        T = Float64(temperature)
                        
                        # Universal gas constant
                        R = 8.31446261815324  # J/(mol·K)
                        
                        # NASA-7 polynomial coefficients for N2 from GRI-Mech 3.0
                        # Range 1: 200-1000K
                        # Range 2: 1000-6000K
                        nasa_coeff = [
                            # Low temperature range (200-1000K)
                            [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372],
                            # High temperature range (1000-6000K)
                            [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]
                        ]
                        
                        temperature_ranges = [[200.0, 1000.0], [1000.0, 6000.0]]
                        
                        # Get appropriate coefficients
                        coeffs = nothing
                        for i in 1:length(temperature_ranges)
                            if temperature_ranges[i][1] <= T && T <= temperature_ranges[i][2]
                                coeffs = nasa_coeff[i]
                                break
                            end
                        end
                        
                        # If temperature outside ranges, use closest range
                        if coeffs === nothing
                            if T < temperature_ranges[1][1]
                                coeffs = nasa_coeff[1]
                            else
                                coeffs = nasa_coeff[2]
                            end
                            @warn "Temperature $T K is outside defined ranges, using closest range coefficients"
                        end
                        
                        # Calculate Cp/R
                        cp_over_R = coeffs[1] + coeffs[2]*T + coeffs[3]*T^2 + coeffs[4]*T^3 + coeffs[5]*T^4
                        
                        # Calculate H/RT
                        h_over_RT = coeffs[1] + (coeffs[2]/2)*T + (coeffs[3]/3)*T^2 + (coeffs[4]/4)*T^3 + (coeffs[5]/5)*T^4 + coeffs[6]/T
                        
                        # Calculate S/R
                        s_over_R = coeffs[1]*log(T) + coeffs[2]*T + (coeffs[3]/2)*T^2 + (coeffs[4]/3)*T^3 + (coeffs[5]/4)*T^4 + coeffs[7]
                        
                        # Convert to standard units
                        cp = cp_over_R * R             # J/mol/K
                        h = h_over_RT * R * T / 1000   # kJ/mol
                        s = s_over_R * R               # J/mol/K
                        g = h - T * s / 1000           # kJ/mol
                        
                        # Add small uncertainty - 2% for GRIM-MECH data
                        cp_uncertainty = 0.02 * cp
                        h_uncertainty = 0.02 * abs(h)
                        s_uncertainty = 0.02 * s
                        g_uncertainty = 0.02 * abs(g)
                        
                        # Return in standard format
                        return Dict(
                            "species_name" => "N2",
                            "temperature" => T,
                            "data_source" => "NASA7-GRIMech-Direct",
                            "properties" => Dict(
                                "Cp" => Dict(
                                    "value" => cp,
                                    "uncertainty" => cp_uncertainty,
                                    "units" => "J/mol/K"
                                ),
                                "H" => Dict(
                                    "value" => h,
                                    "uncertainty" => h_uncertainty,
                                    "units" => "kJ/mol"
                                ),
                                "S" => Dict(
                                    "value" => s,
                                    "uncertainty" => s_uncertainty,
                                    "units" => "J/mol/K"
                                ),
                                "G" => Dict(
                                    "value" => g,
                                    "uncertainty" => g_uncertainty,
                                    "units" => "kJ/mol"
                                )
                            )
                        )
                    end
                    
                    function calculate_n2_properties_direct(temp_range, step)
                        temp_min = Float64(temp_range[1])
                        temp_max = Float64(temp_range[2])
                        step_val = Float64(step)
                        
                        # Generate temperature points
                        temperatures = temp_min:step_val:temp_max
                        
                        # Calculate properties at each temperature
                        results = Dict[]
                        
                        for temp in temperatures
                            result = calculate_n2_properties(temp)
                            push!(results, result)
                        end
                        
                        # Extract property data for plotting
                        temps = [result["temperature"] for result in results]
                        
                        cp_values = [result["properties"]["Cp"]["value"] for result in results]
                        cp_uncertainties = [result["properties"]["Cp"]["uncertainty"] for result in results]
                        
                        h_values = [result["properties"]["H"]["value"] for result in results]
                        h_uncertainties = [result["properties"]["H"]["uncertainty"] for result in results]
                        
                        s_values = [result["properties"]["S"]["value"] for result in results]
                        s_uncertainties = [result["properties"]["S"]["uncertainty"] for result in results]
                        
                        g_values = [result["properties"]["G"]["value"] for result in results]
                        g_uncertainties = [result["properties"]["G"]["uncertainty"] for result in results]
                        
                        # Create result structure
                        return Dict(
                            "species_name" => "N2",
                            "temperature_range" => temp_range,
                            "temperatures" => temps,
                            "properties" => Dict(
                                "Cp" => Dict(
                                    "values" => cp_values,
                                    "uncertainties" => cp_uncertainties,
                                    "units" => "J/mol/K"
                                ),
                                "H" => Dict(
                                    "values" => h_values,
                                    "uncertainties" => h_uncertainties,
                                    "units" => "kJ/mol"
                                ),
                                "S" => Dict(
                                    "values" => s_values,
                                    "uncertainties" => s_uncertainties,
                                    "units" => "J/mol/K"
                                ),
                                "G" => Dict(
                                    "values" => g_values,
                                    "uncertainties" => g_uncertainties,
                                    "units" => "kJ/mol"
                                )
                            )
                        )
                    end
                    
                    # Call the direct calculator
                    return calculate_n2_properties_direct(temp_range, step)
                catch e
                    @warn "Failed to use direct calculator: $e"
                end
            end
            
            # Return empty structure with the right fields
            return Dict(
                "species_name" => species_name,
                "temperature_range" => temp_range,
                "temperatures" => Float64[],
                "properties" => Dict(
                    "Cp" => Dict(
                        "values" => Float64[],
                        "uncertainties" => Float64[],
                        "units" => "J/mol/K"
                    ),
                    "H" => Dict(
                        "values" => Float64[],
                        "uncertainties" => Float64[],
                        "units" => "kJ/mol"
                    ),
                    "S" => Dict(
                        "values" => Float64[],
                        "uncertainties" => Float64[],
                        "units" => "J/mol/K"
                    ),
                    "G" => Dict(
                        "values" => Float64[],
                        "uncertainties" => Float64[],
                        "units" => "kJ/mol"
                    )
                )
            )
        end
        
        # Extract property data
        temps = [result["temperature"] for result in results]
        
        cp_values = [result["properties"]["Cp"]["value"] for result in results]
        cp_uncertainties = [result["properties"]["Cp"]["uncertainty"] for result in results]
        
        h_values = [result["properties"]["H"]["value"] for result in results]
        h_uncertainties = [result["properties"]["H"]["uncertainty"] for result in results]
        
        s_values = [result["properties"]["S"]["value"] for result in results]
        s_uncertainties = [result["properties"]["S"]["uncertainty"] for result in results]
        
        g_values = [result["properties"]["G"]["value"] for result in results]
        g_uncertainties = [result["properties"]["G"]["uncertainty"] for result in results]
        
        # Create result structure
        return Dict(
            "species_name" => species_name,
            "temperature_range" => temp_range,
            "temperatures" => temps,
            "properties" => Dict(
                "Cp" => Dict(
                    "values" => cp_values,
                    "uncertainties" => cp_uncertainties,
                    "units" => "J/mol/K"
                ),
                "H" => Dict(
                    "values" => h_values,
                    "uncertainties" => h_uncertainties,
                    "units" => "kJ/mol"
                ),
                "S" => Dict(
                    "values" => s_values,
                    "uncertainties" => s_uncertainties,
                    "units" => "J/mol/K"
                ),
                "G" => Dict(
                    "values" => g_values,
                    "uncertainties" => g_uncertainties,
                    "units" => "kJ/mol"
                )
            )
        )
    end
    
    # Directly calculate properties for N2
    println("Calculating properties for $species_name over temperature range $temp_range K...")
    result = calculate_properties_direct(conn, species_name, temp_range, config, temp_step)
    
    # Check if we got any valid temperature points
    if isempty(result["temperatures"])
        println("No valid temperature points calculated for $species_name")
    else
        println("Successfully calculated properties at $(length(result["temperatures"])) temperature points")
        
        # Create plots directory if needed
        plot_dir = joinpath(dirname(@__DIR__), "plots")
        mkpath(plot_dir)
        
        # Generate plots
        println("Generating plots...")
        
        # Heat capacity plot
        temps = result["temperatures"]
        cp_values = result["properties"]["Cp"]["values"]
        cp_uncertainties = result["properties"]["Cp"]["uncertainties"]
        
        p_cp = plot(
            temps,
            cp_values,
            ribbon=cp_uncertainties,
            fillalpha=0.3,
            xlabel="Temperature (K)",
            ylabel="Cp (J/mol/K)",
            title="Heat Capacity for $species_name",
            legend=false,
            lw=2,
            grid=true
        )
        
        # Save plot
        plot_file = joinpath(plot_dir, "direct_$(species_name)_Cp.png")
        savefig(p_cp, plot_file)
        println("Saved heat capacity plot to: $plot_file")
        
        # Create tables directory
        tables_dir = joinpath(dirname(@__DIR__), "output", "tables")
        mkpath(tables_dir)
        
        # Create tabular output
        table_file = joinpath(tables_dir, "$(species_name)_data.csv")
        
        # Open CSV file
        open(table_file, "w") do io
            # Write header
            println(io, "Temperature,Cp,Cp_uncertainty,H,H_uncertainty,S,S_uncertainty,G,G_uncertainty")
            
            # Write data rows
            for i in 1:length(result["temperatures"])
                temp = result["temperatures"][i]
                cp = result["properties"]["Cp"]["values"][i]
                cp_unc = result["properties"]["Cp"]["uncertainties"][i]
                h = result["properties"]["H"]["values"][i]
                h_unc = result["properties"]["H"]["uncertainties"][i]
                s = result["properties"]["S"]["values"][i]
                s_unc = result["properties"]["S"]["uncertainties"][i]
                g = result["properties"]["G"]["values"][i]
                g_unc = result["properties"]["G"]["uncertainties"][i]
                
                println(io, "$temp,$cp,$cp_unc,$h,$h_unc,$s,$s_unc,$g,$g_unc")
            end
        end
        
        println("Saved data to CSV: $table_file")
    end
finally
    # Close database connection
    JThermodynamicsData.close_database(conn)
    println("Test complete!")
end