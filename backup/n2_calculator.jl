#!/usr/bin/env julia

"""
N2 Thermodynamic Properties Calculator

This script calculates thermodynamic properties for N2 using NASA-7 polynomial 
coefficients directly, without requiring database access.
"""

# Activate the project to ensure we use the right dependencies
using Pkg
Pkg.activate(@__DIR__)

using Plots
using Printf
using Dates

"""
    calculate_n2_properties(temperature::Real)

Direct calculator for N2 using NASA-7 polynomial coefficients from GRI-Mech 3.0.
"""
function calculate_n2_properties(temperature::Real)
    # Convert temperature to Float64
    T = Float64(temperature)
    
    # Universal gas constant
    R = 8.31446261815324  # J/(mol·K)
    
    # NASA-7 polynomial coefficients for N2 from GRI-Mech 3.0
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

Calculate N2 properties over a temperature range using direct NASA-7 polynomials.
"""
function calculate_n2_properties_range(temp_range::Vector{<:Real}, step::Real=10.0)
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

function main()
    println("N2 Thermodynamic Properties Calculator")
    println("======================================")
    println("Start time: $(now())")
    println()
    
    # Set temperature range and step size
    temp_range = [100.0, 5000.0]
    temp_step = 100.0
    
    println("Calculating properties for N2 over temperature range $(temp_range[1])-$(temp_range[2]) K")
    println("Using temperature step: $temp_step K")
    
    # Calculate properties
    result = calculate_n2_properties_range(temp_range, temp_step)
    
    # Create plots directory if needed
    plot_dir = joinpath(@__DIR__, "plots")
    mkpath(plot_dir)
    
    # Generate plots
    println("\nGenerating plots...")
    
    # Set up plot parameters
    temps = result["temperatures"]
    
    # Heat capacity plot
    p_cp = plot(
        temps,
        result["properties"]["Cp"]["values"],
        ribbon=result["properties"]["Cp"]["uncertainties"],
        fillalpha=0.3,
        xlabel="Temperature (K)",
        ylabel="Cp (J/mol/K)",
        title="Heat Capacity for N2",
        label="NASA-7 Polynomial",
        lw=2,
        grid=true
    )
    
    cp_file = joinpath(plot_dir, "N2_Cp.png")
    savefig(p_cp, cp_file)
    println("  - Heat capacity plot saved to: $cp_file")
    
    # Enthalpy plot
    p_h = plot(
        temps,
        result["properties"]["H"]["values"],
        ribbon=result["properties"]["H"]["uncertainties"],
        fillalpha=0.3,
        xlabel="Temperature (K)",
        ylabel="H (kJ/mol)",
        title="Enthalpy for N2",
        label="NASA-7 Polynomial",
        lw=2,
        grid=true
    )
    
    h_file = joinpath(plot_dir, "N2_H.png")
    savefig(p_h, h_file)
    println("  - Enthalpy plot saved to: $h_file")
    
    # Entropy plot
    p_s = plot(
        temps,
        result["properties"]["S"]["values"],
        ribbon=result["properties"]["S"]["uncertainties"],
        fillalpha=0.3,
        xlabel="Temperature (K)",
        ylabel="S (J/mol/K)",
        title="Entropy for N2",
        label="NASA-7 Polynomial",
        lw=2,
        grid=true
    )
    
    s_file = joinpath(plot_dir, "N2_S.png")
    savefig(p_s, s_file)
    println("  - Entropy plot saved to: $s_file")
    
    # Gibbs energy plot
    p_g = plot(
        temps,
        result["properties"]["G"]["values"],
        ribbon=result["properties"]["G"]["uncertainties"],
        fillalpha=0.3,
        xlabel="Temperature (K)",
        ylabel="G (kJ/mol)",
        title="Gibbs Energy for N2",
        label="NASA-7 Polynomial",
        lw=2,
        grid=true
    )
    
    g_file = joinpath(plot_dir, "N2_G.png")
    savefig(p_g, g_file)
    println("  - Gibbs energy plot saved to: $g_file")
    
    # Combined plot
    p_combined = plot(p_cp, p_h, p_s, p_g, layout=(2,2), size=(1000, 800))
    combined_file = joinpath(plot_dir, "N2_combined.png")
    savefig(p_combined, combined_file)
    println("  - Combined properties plot saved to: $combined_file")
    
    # Create tables directory
    tables_dir = joinpath(@__DIR__, "output", "tables")
    mkpath(tables_dir)
    
    # Create CSV file
    table_file = joinpath(tables_dir, "N2_data.csv")
    
    println("\nCreating tabular output...")
    open(table_file, "w") do io
        # Write header
        println(io, "Temperature,Cp,Cp_uncertainty,H,H_uncertainty,S,S_uncertainty,G,G_uncertainty")
        
        # Write data rows
        for i in 1:length(temps)
            T = temps[i]
            cp = result["properties"]["Cp"]["values"][i]
            cp_unc = result["properties"]["Cp"]["uncertainties"][i]
            h = result["properties"]["H"]["values"][i]
            h_unc = result["properties"]["H"]["uncertainties"][i]
            s = result["properties"]["S"]["values"][i]
            s_unc = result["properties"]["S"]["uncertainties"][i]
            g = result["properties"]["G"]["values"][i]
            g_unc = result["properties"]["G"]["uncertainties"][i]
            
            println(io, "$T,$cp,$cp_unc,$h,$h_unc,$s,$s_unc,$g,$g_unc")
        end
    end
    
    println("  - Data saved to: $table_file")
    
    println("\nCalculation and visualization complete!")
    println("End time: $(now())")
end

# Run the main function
main()