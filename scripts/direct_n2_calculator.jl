#!/usr/bin/env julia

# Direct N2 thermodynamic property calculator
# This script provides a simplified implementation for calculating
# thermodynamic properties of N2 using the NASA-7 polynomial

"""
    calculate_n2_properties(temperature::Real)

Calculate thermodynamic properties of N2 at a specified temperature using 
NASA-7 polynomial coefficients from GRI-Mech 3.0.
"""
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

"""
    calculate_n2_properties_range(temp_range::Vector{<:Real}, step::Real=10.0)

Calculate thermodynamic properties of N2 over a temperature range.
"""
function calculate_n2_properties_range(temp_range::Vector{<:Real}, step::Real=10.0)
    temp_min = Float64(temp_range[1])
    temp_max = Float64(temp_range[2])
    temp_step = Float64(step)
    
    # Generate temperature points
    temperatures = temp_min:temp_step:temp_max
    
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

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    println("Direct N2 Property Calculator")
    println("============================")
    
    # Calculate at standard temperature
    result = calculate_n2_properties(298.15)
    
    println("\nProperties for N2 at 298.15 K:")
    for (prop, data) in result["properties"]
        println("$prop: $(data["value"]) ± $(data["uncertainty"]) $(data["units"])")
    end
    
    # Calculate over a temperature range
    temp_range = [200.0, 3000.0]
    step = 100.0
    range_result = calculate_n2_properties_range(temp_range, step)
    
    println("\nCalculated properties over temperature range $(temp_range[1])-$(temp_range[2]) K")
    println("Temperature points: $(length(range_result["temperatures"]))")
    
    # Try to create plots if Plots package is available
    try
        using Plots
        
        # Create plots
        plot_dir = joinpath(dirname(@__DIR__), "plots")
        mkpath(plot_dir)
        
        # Heat capacity plot
        temps = range_result["temperatures"]
        cp_values = range_result["properties"]["Cp"]["values"]
        cp_uncertainties = range_result["properties"]["Cp"]["uncertainties"]
        
        p_cp = plot(
            temps,
            cp_values,
            ribbon=cp_uncertainties,
            fillalpha=0.3,
            xlabel="Temperature (K)",
            ylabel="Cp (J/mol/K)",
            title="Heat Capacity for N2",
            legend=false,
            lw=2,
            grid=true
        )
        
        # Save plot
        plot_file = joinpath(plot_dir, "direct_N2_Cp.png")
        savefig(p_cp, plot_file)
        println("\nSaved heat capacity plot to: $plot_file")
        
        # Entropy plot
        s_values = range_result["properties"]["S"]["values"]
        s_uncertainties = range_result["properties"]["S"]["uncertainties"]
        
        p_s = plot(
            temps,
            s_values,
            ribbon=s_uncertainties,
            fillalpha=0.3,
            xlabel="Temperature (K)",
            ylabel="S (J/mol/K)",
            title="Entropy for N2",
            legend=false,
            lw=2,
            grid=true
        )
        
        # Save plot
        plot_file = joinpath(plot_dir, "direct_N2_S.png")
        savefig(p_s, plot_file)
        println("Saved entropy plot to: $plot_file")
    catch e
        println("\nPlotting skipped: Plots package not available or error occurred: $e")
        println("You can install Plots with: using Pkg; Pkg.add(\"Plots\")")
    end
    
    # Save data as CSV for verification
    outdir = joinpath(dirname(@__DIR__), "output", "tables")
    mkpath(outdir)
    
    csv_file = joinpath(outdir, "direct_N2_data.csv")
    open(csv_file, "w") do io
        # Write header
        println(io, "Temperature,Cp,Cp_uncertainty,H,H_uncertainty,S,S_uncertainty,G,G_uncertainty")
        
        # Write data rows
        temps = range_result["temperatures"]
        n = length(temps)
        
        for i in 1:n
            temp = temps[i]
            cp = range_result["properties"]["Cp"]["values"][i]
            cp_unc = range_result["properties"]["Cp"]["uncertainties"][i]
            h = range_result["properties"]["H"]["values"][i]
            h_unc = range_result["properties"]["H"]["uncertainties"][i]
            s = range_result["properties"]["S"]["values"][i]
            s_unc = range_result["properties"]["S"]["uncertainties"][i]
            g = range_result["properties"]["G"]["values"][i]
            g_unc = range_result["properties"]["G"]["uncertainties"][i]
            
            println(io, "$temp,$cp,$cp_unc,$h,$h_unc,$s,$s_unc,$g,$g_unc")
        end
    end
    
    println("\nSaved data to CSV: $csv_file")
    println("\nThis direct calculator can be used to debug and verify the main package functionality.")
end