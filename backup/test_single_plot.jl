#\!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using Printf
using Plots
using JSON
using Colors
using Dates

"""
    plot_single_species(species_name)

Plot a single species with detailed debug output to 
verify what sources are being included in the plot.
"""
function plot_single_species(species_name)
    println("\n=== DETAILED PLOT GENERATION FOR $(species_name) ===\n")
    
    # Temperature range
    temp_range = [200.0, 6000.0]
    
    # Get data sources with extra debug output
    println("Loading species data...")
    species_data = JThermodynamicsData.load_species_data(species_name)
    
    if \!haskey(species_data, "sources") || isempty(species_data["sources"])
        println("  ERROR: No data sources available for $(species_name)")
        return false
    end
    
    println("\nDumping full raw source data:")
    for (source_name, source_info) in species_data["sources"]
        priority = get(source_info, "priority", 0)
        reliability = get(source_info, "reliability_score", 0.0)
        println("  * Source: $(source_name) (priority: $(priority), reliability: $(reliability))")
    end
    
    println("\nGetting sources sorted by priority...")
    all_sources = [(name, get(info, "priority", 0), get(info, "reliability_score", 0.0)) 
              for (name, info) in species_data["sources"]]
    sort\!(all_sources, by=s->s[2], rev=true)
    
    println("All sources (sorted by priority):")
    for (src_name, priority, reliability) in all_sources
        println("  - $(src_name) (priority: $(priority), reliability: $(reliability))")
    end
    
    println("\nFiltering sources...")
    # Explicitly filter out all theoretical sources except the best one
    theoretical_sources = filter(s -> startswith(lowercase(s[1]), "theoretical") || lowercase(s[1]) == "theoretical", all_sources)
    experimental_sources = filter(s -> \!startswith(lowercase(s[1]), "theoretical") && lowercase(s[1]) \!= "theoretical", all_sources)
    
    println("Theoretical sources found:")
    for source in theoretical_sources
        println("  * $(source[1]) (priority: $(source[2]), reliability: $(source[3]))")
    end
    
    println("\nExperimental sources found:")
    for source in experimental_sources
        println("  * $(source[1]) (priority: $(source[2]), reliability: $(source[3]))")
    end
    
    # Keep only the best theoretical source
    best_theoretical = isempty(theoretical_sources) ? [] : [theoretical_sources[1]]
    
    println("\nBest theoretical source:")
    for source in best_theoretical
        println("  * $(source[1]) (priority: $(source[2]), reliability: $(source[3]))")
    end
    
    # Replace the sources list with filtered version
    filtered_sources = vcat(best_theoretical, experimental_sources)
    sort\!(filtered_sources, by=s->s[2], rev=true)
    
    println("\nFINAL FILTERED SOURCES (these will be in the plot):")
    for (src_name, priority, reliability) in filtered_sources
        println("  → $(src_name) (priority: $(priority), reliability: $(reliability))")
    end
    
    # Database path
    db_path = joinpath(@__DIR__, "data", "thermodynamics.duckdb")
    
    # Load best source data (will be displayed prominently)
    best_poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, Dict("method" => "hierarchical"))
    best_source = best_poly.source
    
    println("\nSelected source: $(best_source)")
    
    # Map sources to their standardized display name
    source_display_names = Dict()
    for (src_name, _, _) in filtered_sources
        # All theoretical sources get one standard name for consistency
        if startswith(lowercase(src_name), "theoretical") || lowercase(src_name) == "theoretical"
            source_display_names[src_name] = "theoretical"
        else
            # Non-theoretical sources keep their original name
            source_display_names[src_name] = src_name
        end
    end
    
    println("\nDisplay names for plot legend:")
    for (orig_name, display_name) in source_display_names
        println("  * $(orig_name) → '$(display_name)'")
    end
    
    # Generate temperatures on logarithmic scale
    temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=100)
    
    # Create a 2x2 plot layout for the 4 thermodynamic properties with room for legends
    plt = plot(layout=(2,2), size=(1200, 900), dpi=300, margin=10Plots.mm, right_margin=20Plots.mm)
    
    # Properties to plot
    properties = [
        (1, "cp", "Cp (J/mol/K)"),
        (2, "h", "H (kJ/mol)"),
        (3, "s", "S (J/mol/K)"),
        (4, "g", "G (kJ/mol)")
    ]
    
    # Assign colors and line styles
    source_colors = Dict()
    line_styles = Dict()
    
    for (src_name, priority, _) in filtered_sources
        if src_name == best_source
            source_colors[src_name] = :blue
            line_styles[src_name] = :solid
        else
            # All other sources get a gray color based on their priority
            priority_level = max(1, min(priority, 13))
            color_value = 0.7 - (priority_level / 20.0)  # Higher priority = darker color
            source_colors[src_name] = RGB(color_value, color_value, color_value)
            line_styles[src_name] = :dash
        end
    end
    
    # Plot each property for each source in the FILTERED list
    println("\nPlotting each source:")
    for (source_name, priority, reliability) in filtered_sources
        println("  * Creating plot for source: $(source_name)")
        
        try
            # Load data for this source
            poly = JThermodynamicsData.load_polynomial_data(
                db_path, species_name, temp_range, Dict("source" => source_name)
            )
            
            # Set line properties
            line_color = source_colors[source_name]
            line_style = line_styles[source_name]
            line_width = source_name == best_source ? 2.5 : 1.5
            alpha = source_name == best_source ? 1.0 : 0.7
            
            # Use standardized display name for the legend
            display_name = source_display_names[source_name]
            println("    - Using legend name: '$(display_name)'")
            
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
                    
                    plot\!(
                        plt, 
                        temps, 
                        values, 
                        subplot=idx,
                        label=display_name,  # Use standardized display name
                        xlabel="Temperature (K)",
                        ylabel=ylabel,
                        title="$(uppercase(property)) for $(species_name)",
                        xscale=:log10,
                        xticks=([200, 1000, 5000], ["200", "1000", "5000"]),
                        line=(line_color, line_width, line_style),
                        ribbon=(values .- lower, upper .- values),
                        fillalpha=0.2,
                        legend=:outerright
                    )
                else
                    # Other sources - no uncertainty ribbon
                    plot\!(
                        plt, 
                        temps, 
                        values, 
                        subplot=idx,
                        label=display_name,  # Use standardized display name
                        line=(line_color, line_width, line_style),
                        alpha=alpha,
                        legend=:outerright
                    )
                end
            end
        catch e
            println("  Error plotting data for source $(source_name): $(e)")
        end
    end
    
    # Add a main title
    plot\!(plt, title="Thermodynamic Properties for $(species_name)", 
          titleloc=:center, titlefontsize=12)
    
    # Save the plot
    plots_dir = joinpath(@__DIR__, "plots", "debug")
    mkpath(plots_dir)
    plot_file = joinpath(plots_dir, "$(species_name)_debug.png")
    savefig(plt, plot_file)
    
    println("\nDebug plot saved to: $(plot_file)")
    return true
end

# Run the function on a specific problematic species
plot_single_species("OH-")
plot_single_species("Zn")
