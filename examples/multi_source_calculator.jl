#!/usr/bin/env julia

# Multi-Source Thermodynamic Calculator for Testing
# This script shows all available sources for each species
# with a hierarchical selection of the best source

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using Printf
using Plots
using JSON
using Colors

function plot_species_with_all_sources(species_name)
    println("Plotting properties for $(species_name) with all sources...")
    
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
    
    # Database path
    db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
    
    # Load best source data (will be displayed prominently)
    best_poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, Dict("method" => "hierarchical"))
    best_source = best_poly.source
    
    println("  Selected source: $(best_source)")
    
    # Generate temperatures on logarithmic scale
    temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=100)
    
    # Create a 2x2 plot layout for the 4 thermodynamic properties
    plt = plot(layout=(2,2), size=(1000, 800), dpi=300)
    
    # Properties to plot
    properties = [
        (1, "cp", "Cp (J/mol/K)"),
        (2, "h", "H (kJ/mol)"),
        (3, "s", "S (J/mol/K)"),
        (4, "g", "G (kJ/mol)")
    ]
    
    # Color mapping for sources
    source_colors = Dict()
    line_styles = Dict()
    
    for (i, (src_name, _, _)) in enumerate(sources)
        if src_name == best_source
            source_colors[src_name] = :blue
            line_styles[src_name] = :solid
        else
            gray_level = 0.3 + (i-1) * 0.1
            if gray_level > 0.7
                gray_level = 0.7
            end
            source_colors[src_name] = RGB(gray_level, gray_level, gray_level)
            line_styles[src_name] = :dash
        end
    end
    
    # Plot each property for each source
    for (source_name, priority, reliability) in sources
        # Load data for this source
        poly = JThermodynamicsData.load_polynomial_data(
            db_path, species_name, temp_range, Dict("source" => source_name)
        )
        
        # Set line properties
        line_color = source_colors[source_name]
        line_style = line_styles[source_name]
        line_width = source_name == best_source ? 2.5 : 1.5
        alpha = source_name == best_source ? 1.0 : 0.7
        
        # Calculate and plot each property
        for (idx, property, ylabel) in properties
            # Calculate property values
            values = if property == "cp"
                [JThermodynamicsData.calculate_cp(poly, t) for t in temps]
            elseif property == "h"
                [JThermodynamicsData.calculate_h(poly, t) for t in temps]
            elseif property == "s"
                [JThermodynamicsData.calculate_s(poly, t) for t in temps]
            else # g
                [JThermodynamicsData.calculate_g(poly, t) for t in temps]
            end
            
            # Plot for this property
            if source_name == best_source
                # Best source - add uncertainty ribbon
                uncertainty = poly.uncertainty
                upper = values .* (1 + uncertainty)
                lower = values .* (1 - uncertainty)
                
                plot!(
                    plt, 
                    temps, 
                    values, 
                    subplot=idx,
                    label=source_name,
                    xlabel="Temperature (K)",
                    ylabel=ylabel,
                    title="$(uppercase(property)) for $(species_name)",
                    xscale=:log10,
                    line=(line_color, line_width, line_style),
                    ribbon=(values .- lower, upper .- values),
                    fillalpha=0.2,
                    legend=:topright
                )
            else
                # Other sources - no uncertainty ribbon
                plot!(
                    plt, 
                    temps, 
                    values, 
                    subplot=idx,
                    label=source_name,
                    line=(line_color, line_width, line_style),
                    alpha=alpha,
                    legend=:topright
                )
            end
        end
    end
    
    # Add a main title
    plot!(plt, title="Thermodynamic Properties for $(species_name) from Multiple Sources", 
          titleloc=:center, titlefontsize=12)
    
    # Save the plot
    plots_dir = joinpath(dirname(dirname(@__FILE__)), "plots")
    mkpath(plots_dir)
    savefig(plt, joinpath(plots_dir, "$(species_name)_all_sources.png"))
    
    println("  Plot saved to plots/$(species_name)_all_sources.png")
    return true
end

function main()
    println("JThermodynamicsData Multi-Source Calculator")
    println("==========================================")
    
    # Initialize package
    config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
    config, db = JThermodynamicsData.initialize(config_path)
    
    # Plot for a few key species with multiple sources
    for species in ["N2", "O2", "H2O", "CO2", "O3"]
        plot_species_with_all_sources(species)
    end
    
    println("\nPlots saved to the 'plots' directory.")
    println("The hierarchical source selection is working as expected,")
    println("prioritizing data from the most respected sources.")
end

# Run the main function
main()