#!/usr/bin/env julia

# Simple Thermodynamic Calculator for Testing
# This is a simplified version to diagnose plot display issues

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using Printf
using Plots
using JSON

function plot_species_properties(species_name)
    println("Plotting properties for $(species_name)...")
    
    # Temperature range
    temp_range = [200.0, 6000.0]
    
    # Get data sources
    species_data = JThermodynamicsData.load_species_data(species_name)
    
    if !haskey(species_data, "sources") || isempty(species_data["sources"])
        println("  No data sources available for $(species_name)")
        return false
    end
    
    # Get sources and sort by priority
    sources = [(name, get(info, "priority", 0), get(info, "reliability_score", 0.0)) 
              for (name, info) in species_data["sources"]]
    sort!(sources, by=s->s[2], rev=true)
    
    # Print available sources
    println("  Available sources:")
    for (src_name, priority, reliability) in sources
        println("  - $(src_name) (priority: $(priority), reliability: $(reliability))")
    end
    
    # Load best source data
    db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
    best_poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, Dict("method" => "hierarchical"))
    
    println("  Selected source: $(best_poly.source)")
    
    # Generate temperatures on logarithmic scale
    temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=100)
    
    # Calculate properties
    cp_values = [JThermodynamicsData.calculate_cp(best_poly, t) for t in temps]
    h_values = [JThermodynamicsData.calculate_h(best_poly, t) for t in temps]
    s_values = [JThermodynamicsData.calculate_s(best_poly, t) for t in temps]
    g_values = [JThermodynamicsData.calculate_g(best_poly, t) for t in temps]
    
    # Create a simple plot
    plt = plot(
        temps, cp_values, 
        title="Heat Capacity of $(species_name)",
        xlabel="Temperature (K)", 
        ylabel="Cp (J/mol/K)",
        label="$(best_poly.source)",
        xscale=:log10,
        lw=2
    )
    
    # Save the plot
    plots_dir = joinpath(dirname(dirname(@__FILE__)), "plots")
    mkpath(plots_dir)
    savefig(plt, joinpath(plots_dir, "$(species_name)_cp_simple.png"))
    
    println("  Plot saved to plots/$(species_name)_cp_simple.png")
    return true
end

function main()
    println("JThermodynamicsData Simple Calculator")
    println("=====================================")
    
    # Initialize package
    config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
    config, db = JThermodynamicsData.initialize(config_path)
    
    # Plot for a few key species
    for species in ["N2", "O2", "H2O", "CO2", "O3"]
        plot_species_properties(species)
    end
    
    println("\nPlots saved to the 'plots' directory.")
end

# Run the main function
main()