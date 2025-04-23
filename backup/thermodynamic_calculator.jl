#!/usr/bin/env julia

"""
Hierarchical Thermodynamic Properties Calculator

This script uses the JThermodynamicsData package to calculate thermodynamic properties
for chemical species using a hierarchical approach. It dynamically fetches data from
original sources following a hierarchical structure defined in the configuration.

Usage:
    julia thermodynamic_calculator.jl [species] [temperature]

Examples:
    julia thermodynamic_calculator.jl N2 298.15
    julia thermodynamic_calculator.jl H2O 300-1500
"""

# Activate the project to ensure we use the right dependencies
using Pkg
Pkg.activate(@__DIR__)

# Import the JThermodynamicsData package
using JThermodynamicsData
using Plots

"""
    main()

Main entry point for the thermodynamic calculator script.
"""
function main()
    if length(ARGS) < 1
        println("Usage: julia thermodynamic_calculator.jl [species] [temperature]")
        println("Examples:")
        println("  julia thermodynamic_calculator.jl N2 298.15")
        println("  julia thermodynamic_calculator.jl H2O 300-1500")
        return
    end
    
    # Get species and temperature from command line arguments
    species_name = ARGS[1]
    
    temp_range = [200.0, 6000.0]  # Default temperature range
    single_temp = nothing
    
    if length(ARGS) >= 2
        temp_arg = ARGS[2]
        if occursin("-", temp_arg)
            # Range format: "min-max"
            temp_parts = split(temp_arg, "-")
            if length(temp_parts) == 2
                temp_min = parse(Float64, temp_parts[1])
                temp_max = parse(Float64, temp_parts[2])
                temp_range = [temp_min, temp_max]
            end
        else
            # Single temperature - we'll use a small range around it for calculation
            # but display just the single point
            temperature = parse(Float64, temp_arg)
            single_temp = temperature
            temp_range = [max(temperature - 10, 200.0), min(temperature + 10, 6000.0)]
        end
    end
    
    # Make sure the data directory exists
    mkpath(joinpath(@__DIR__, "data"))
    
    # Initialize the package with default configuration
    config, db = initialize()
    
    println("Calculating thermodynamic properties for $species_name over temperature range $(temp_range[1])-$(temp_range[2]) K")
    
    try
        # Use the hierarchical calculator to get thermodynamic data
        result = calculate_properties(db, species_name, temp_range, Dict("method" => "hierarchical"))
        
        # Display the results
        if single_temp === nothing
            # Plot the results for a temperature range
            plots = plot_properties(result)
            
            # Create output directory if it doesn't exist
            mkpath(joinpath(@__DIR__, "output"))
            
            for (prop, plt) in plots
                display(plt)
                savefig(plt, joinpath(@__DIR__, "output", "$(species_name)_$(prop).png"))
            end
        else
            # Display numerical results for a single temperature
            println("\nThermodynamic properties at $single_temp K:")
            println("============================================")
            
            # Check if we have any data 
            if !haskey(result, "temperatures") || isempty(result["temperatures"])
                println("No data found for species $species_name")
                return
            end
            
            # Get the result for the point closest to our requested temperature
            temps = result["temperatures"]
            idx = argmin(abs.(temps .- single_temp))
            closest_temp = temps[idx]
            
            props = ["Cp", "H", "S", "G"]
            for prop in props
                values = result[lowercase(prop)]
                uncertainties = result["$(lowercase(prop))_uncertainty"]
                
                value = values[idx]
                uncertainty = uncertainties[idx]
                
                units = prop == "Cp" ? "J/mol·K" : prop == "H" ? "kJ/mol" : prop == "S" ? "J/mol·K" : "kJ/mol"
                
                println("$prop: $value ± $uncertainty $units")
            end
        end
        
        # Print the sources used (highest priority source first)
        if haskey(result, "sources_used") && !isempty(result["sources_used"])
            highest_source = result["highest_priority_source"]
            println("\nData sources used (in priority order):")
            println("======================================")
            for source in reverse(sort(result["sources_used"], by=x->x["priority"]))
                mark = source["source_name"] == highest_source ? "✓" : " "
                println("$mark $(source["source_name"]) (priority: $(source["priority"]), reliability: $(source["reliability_score"]))")
            end
        else
            println("\nNo data sources were found for $species_name")
        end
    catch e
        println("\nError calculating properties for $species_name: $e")
        # Add suggestions for next steps
        println("\nSuggestions:")
        println("1. Check if the species name is correct")
        println("2. Check if the temperature is within supported range (200-6000K)")
        println("3. Check if data sources are properly configured")
        println("4. Try initializing the database with: ./scripts/initialize_database.jl")
    end
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end