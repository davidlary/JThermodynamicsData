#!/usr/bin/env julia

# Direct Thermodynamic Calculator Demo
# This script demonstrates using the proper JThermodynamicsData calculator
# with dynamic data fetching from original sources
# It creates plots with logarithmic temperature scales and uncertainty shading

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using Printf
using Plots
using Colors
using YAML
using DataFrames
using DuckDB
using CSV
using JSON

"""
    get_all_species(db)

Get list of all species from the database.
"""
function get_all_species(db)
    # First try using JSON storage
    json_species = JThermodynamicsData.list_available_species()
    if !isempty(json_species)
        return json_species
    end
    
    # Fallback to database
    conn = JThermodynamicsData.get_connection(db)
    result = DBInterface.execute(conn, "SELECT name FROM species ORDER BY name")
    df = DataFrame(result)
    
    if size(df, 1) == 0
        # Fallback to hardcoded list if database is empty
        return [
            "N2", "O2", "H2O", "CO2", "CO", "CH4", "H2", "NO", "NO2", "Ar", "He", "NH3", "OH",
            "N", "O", "H", "C", "S", "Xe", "Ne", "SO", "SO2", "SO3", "CN", "HCN", "NCO", "NCN",
            "HNO", "HO2", "H2O2", "HONO", "HNO3", "N2O", "N2H4", "CH2", "CH3", "CHO", "CH2O",
            "CHOOH", "HCNO", "N+", "O+", "H+", "C+", "S+", "Xe+", "He+", "CO+", "NO+", "CN+",
            "N2+", "O2+", "H2+", "CO2+", "NO2+", "N2O+", "H2O+", "H3O+", "CHO+", "HNO+",
            "HNO2+", "HNO3+", "O3+", "N-", "O-", "H-", "C-", "CN-", "CO-", "NO-", "O2-",
            "O3", "OH+", "OH-", "HO2-", "NO2-", "NO3-", "e-"
        ]
    end
    
    return df.name
end

"""
    export_property_data(db, species, temp_range=[200.0, 6000.0])

Export thermodynamic property data for a species to CSV files.
"""
function export_property_data(db, species, temp_range=[200.0, 6000.0])
    try
        # Load data using hierarchical selection
        poly = JThermodynamicsData.load_polynomial_data(db, species, temp_range, Dict("method" => "hierarchical"))
        
        # Generate temperature points (using logarithmic scale)
        temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=200)
        
        # Calculate all properties
        cp_values = [JThermodynamicsData.calculate_cp(poly, t) for t in temps]
        h_values = [JThermodynamicsData.calculate_h(poly, t) for t in temps]
        s_values = [JThermodynamicsData.calculate_s(poly, t) for t in temps]
        g_values = [JThermodynamicsData.calculate_g(poly, t) for t in temps]
        
        # Create a DataFrame
        df = DataFrame(
            Temperature = temps,
            Cp = cp_values,
            H = h_values,
            S = s_values,
            G = g_values,
            Source = fill(poly.source, length(temps)),
            Uncertainty = fill(poly.uncertainty, length(temps))
        )
        
        # Create output directory
        tables_dir = joinpath(dirname(dirname(@__FILE__)), "output", "tables")
        mkpath(tables_dir)
        
        # Export to CSV
        CSV.write(joinpath(tables_dir, "$(species)_hierarchical.csv"), df)
        
        # Also export the polynomial coefficients
        coeff_data = Dict(
            "species" => species,
            "source" => poly.source,
            "reliability_score" => poly.reliability_score,
            "polynomial_type" => poly.polynomial_type,
            "temperature_ranges" => poly.temperature_ranges,
            "coefficients" => poly.coefficients,
            "uncertainty" => poly.uncertainty
        )
        
        open(joinpath(tables_dir, "$(species)_nasa7.txt"), "w") do f
            JSON.print(f, coeff_data, 4)  # Pretty print with 4-space indent
        end
        
        # Create documentation about the source selection
        docs_dir = joinpath(dirname(dirname(@__FILE__)), "output", "docs")
        mkpath(docs_dir)
        
        # JSON report
        open(joinpath(docs_dir, "$(species)_sources.json"), "w") do f
            JSON.print(f, Dict(
                "species" => species,
                "selected_source" => poly.source,
                "reliability_score" => poly.reliability_score,
                "uncertainty" => poly.uncertainty,
                "temperature_range" => temp_range,
                "polynomial_type" => poly.polynomial_type
            ), 4)
        end
        
        # Markdown report
        open(joinpath(docs_dir, "$(species)_sources.md"), "w") do f
            write(f, "# Thermodynamic Data Source for $(species)\n\n")
            write(f, "## Selected Source\n")
            write(f, "- **Source:** $(poly.source)\n")
            write(f, "- **Reliability score:** $(poly.reliability_score)\n")
            write(f, "- **Uncertainty:** $(poly.uncertainty * 100)%\n")
            write(f, "- **Temperature range:** $(temp_range[1]) - $(temp_range[2]) K\n")
            write(f, "- **Polynomial type:** $(poly.polynomial_type)\n\n")
            write(f, "## Coefficients\n")
            write(f, "### Low temperature range ($(poly.temperature_ranges[1][1]) - $(poly.temperature_ranges[1][2]) K)\n")
            write(f, "```\n$(join([@sprintf("%.8e", c) for c in poly.coefficients[1]], ", "))\n```\n\n")
            write(f, "### High temperature range ($(poly.temperature_ranges[2][1]) - $(poly.temperature_ranges[2][2]) K)\n")
            write(f, "```\n$(join([@sprintf("%.8e", c) for c in poly.coefficients[2]], ", "))\n```\n")
            
            # Add information about all available sources
            if haskey(JThermodynamicsData.load_species_data(species), "sources")
                sources = JThermodynamicsData.list_available_sources(species)
                if !isempty(sources)
                    write(f, "\n## All Available Sources\n")
                    for (src_name, priority, reliability) in sources
                        write(f, "- $(src_name) (priority: $(priority), reliability: $(reliability))\n")
                    end
                end
            end
        end
        
        return true
    catch e
        println("  Error exporting data for $(species): $(e)")
        return false
    end
end

"""
    plot_all_properties(db, species, temp_range=[200.0, 6000.0])

Create a plot with all 4 thermodynamic properties for a species.
Uses logarithmic temperature scales and uncertainty shading.
Shows data from all available sources with the highest priority one highlighted.
"""
function plot_all_properties(db, species, temp_range=[200.0, 6000.0])
    # Set up a 2x2 subplot layout with larger margins for better appearance
    plt = plot(
        layout=(2,2),
        size=(1200, 900),
        dpi=300,
        margin=8Plots.mm,
        fontfamily="Computer Modern",
        titlefontsize=10,
        labelfontsize=9
    )
    
    # Define properties to plot and their positions
    properties = [
        (1, "cp", "Cp (J/mol/K)"),
        (2, "h", "H (kJ/mol)"),
        (3, "s", "S (J/mol/K)"),
        (4, "g", "G (kJ/mol)")
    ]
    
    try
        # Get all available sources for this species from JSON storage
        species_data = JThermodynamicsData.load_species_data(species)
        
        if !haskey(species_data, "sources") || isempty(species_data["sources"])
            println("  No data sources available for $(species)")
            return false
        end
        
        # Get sources and sort by priority (highest first)
        sources = [(name, get(info, "priority", 0), get(info, "reliability_score", 0.0)) 
                  for (name, info) in species_data["sources"]]
        sort!(sources, by=s->s[2], rev=true)
        
        # Generate temperature points on logarithmic scale
        temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=300)
        
        # Load the best source data (hierarchical selection)
        best_poly = JThermodynamicsData.load_polynomial_data(db, species, temp_range, Dict("method" => "hierarchical"))
        best_source = best_poly.source
        
        # Export the data (saves both CSV and polynomial coefficients)
        export_property_data(db, species, temp_range)
        
        # Add a main title to the plot
        plot!(plt, title="Thermodynamic Properties for $(species) - Multiple Sources", 
              titlefontsize=14)
        
        # Create a color palette for different sources
        source_colors = Dict()
        line_styles = Dict()
        
        # Assign colors and line styles
        if length(sources) == 1
            # Only one source - use blue
            source_colors[sources[1][1]] = :blue
            line_styles[sources[1][1]] = :solid
        else
            # Multiple sources - use a color gradient
            for (i, (src_name, _, _)) in enumerate(sources)
                if src_name == best_source
                    # Best source - solid blue
                    source_colors[src_name] = :blue
                    line_styles[src_name] = :solid
                else
                    # Other sources - graded from dark gray to light gray
                    gray_level = 0.3 + (i-1) * 0.1
                    if gray_level > 0.7
                        gray_level = 0.7
                    end
                    source_colors[src_name] = RGB(gray_level, gray_level, gray_level)
                    line_styles[src_name] = :dash
                end
            end
        end
        
        # Create legend items for each subplot
        legend_items = Dict()
        for idx in 1:4
            legend_items[idx] = []
        end
        
        # Plot each property for each source
        for (source_name, priority, reliability) in sources
            try
                # Load polynomial data for this source
                poly = JThermodynamicsData.load_polynomial_data(db, species, temp_range, Dict("source" => source_name))
                
                # Get line properties
                line_color = source_colors[source_name]
                line_style = line_styles[source_name]
                line_width = source_name == best_source ? 2.5 : 1.5
                alpha = source_name == best_source ? 1.0 : 0.7
                
                # For each property type
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
                    
                    # Create plot for this property
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
                            title="$(uppercase(property))",
                            xscale=:log10,
                            line=(line_color, line_width, line_style),
                            ribbon=(values .- lower, upper .- values),
                            fillalpha=0.2,
                            grid=true,
                            minorgrid=true,
                            guidefontsize=9
                        )
                        
                        # Add annotation for uncertainty
                        if idx == 1
                            annotate!(
                                plt,
                                [(temp_range[1]*1.1, maximum(upper)*0.95, 
                                text("Uncertainty: ±$(round(uncertainty*100, digits=1))%", 8, :left, :top, :black))]
                            )
                        end
                    else
                        # Other sources - no uncertainty ribbon
                        plot!(
                            plt, 
                            temps, 
                            values, 
                            subplot=idx,
                            label=source_name,
                            line=(line_color, line_width, line_style),
                            alpha=alpha
                        )
                    end
                    
                    # Add to legend items
                    push!(legend_items[idx], (source_name, line_color, line_width, line_style, alpha))
                end
            catch e
                println("  Error plotting data for $(species) from source $(source_name): $(e)")
            end
        end
        
        # Add legends to the plots
        for idx in 1:4
            if !isempty(legend_items[idx])
                # Sort by priority (best source first)
                sort!(legend_items[idx], by=x->x[1] == best_source ? 0 : 1)
                
                # Add source legend to first subplot only
                if idx == 1
                    # Create legend text
                    legend_text = "Data Sources:"
                    for (src_name, _, _, _, _) in legend_items[idx]
                        suffix = src_name == best_source ? " (best)" : ""
                        legend_text *= "\n$(src_name)$(suffix)"
                    end
                    
                    # Position in top right corner
                    y_pos = 0
                    try
                        y_pos = maximum([p[:y][end] for p in plt[idx]])
                    catch
                        # Fallback if there's an error accessing the subplot
                        y_pos = temp_range[2] * 0.1
                    end
                    
                    if y_pos == 0 || isnan(y_pos)
                        y_pos = temp_range[2] * 0.1  # Another fallback
                    end
                    
                    # Add annotation
                    annotate!(
                        plt,
                        [(temp_range[2]*0.95, y_pos*0.9, 
                        text(legend_text, 8, :right, :top, :black))],
                        subplot=idx
                    )
                end
            end
        end
        
        # Save the plot
        plots_dir = joinpath(dirname(dirname(@__FILE__)), "plots")
        mkpath(plots_dir)
        savefig(plt, joinpath(plots_dir, "$(species)_all_properties.png"))
        
        return true
    catch e
        println("  Error processing $(species): $(e)")
        return false
    end
end

"""
    create_summary_report(db, all_species, success_count)

Create a summary report of all species processed.
"""
function create_summary_report(db, all_species, success_count, success_species)
    tables_dir = joinpath(dirname(dirname(@__FILE__)), "output", "tables")
    mkpath(tables_dir)
    
    # Create a summary DataFrame
    species_df = DataFrame(
        Species = all_species,
        Processed = [s in success_species for s in all_species]
    )
    
    # Add additional information if available
    sources = String[]
    uncertainties = Float64[]
    
    for species in all_species
        if species in success_species
            # Try to read the JSON source file
            json_path = joinpath(dirname(dirname(@__FILE__)), "output", "docs", "$(species)_sources.json")
            if isfile(json_path)
                data = JSON.parsefile(json_path)
                push!(sources, get(data, "selected_source", "unknown"))
                push!(uncertainties, get(data, "uncertainty", 0.0))
            else
                push!(sources, "unknown")
                push!(uncertainties, 0.0)
            end
        else
            push!(sources, "failed")
            push!(uncertainties, NaN)
        end
    end
    
    # Add to dataframe
    species_df.Source = sources
    species_df.Uncertainty = uncertainties
    
    # Write summary to CSV
    CSV.write(joinpath(tables_dir, "all_species_summary.csv"), species_df)
    
    return species_df
end

# Main function for the demo
function main()
    println("JThermodynamicsData Direct Calculator Demo")
    println("===========================================")
    
    # Initialize package with configuration
    config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
    config, db = JThermodynamicsData.initialize(config_path)
    
    # Test hierarchical selection for N2
    println("\nTesting hierarchical selection for N2:")
    temp = 298.15
    
    try
        # Use hierarchical selection
        poly = JThermodynamicsData.load_polynomial_data(db, "N2", [temp-10, temp+10], Dict("method" => "hierarchical"))
        
        # Calculate properties
        cp = JThermodynamicsData.calculate_cp(poly, temp)
        h = JThermodynamicsData.calculate_h(poly, temp)
        s = JThermodynamicsData.calculate_s(poly, temp)
        g = JThermodynamicsData.calculate_g(poly, temp)
        
        println("  Selected source: $(poly.source)")
        println("  Properties for N2 at $(temp) K:")
        println("  Cp = $(@sprintf("%.2f", cp)) J/mol/K")
        println("  H  = $(@sprintf("%.2f", h)) kJ/mol")
        println("  S  = $(@sprintf("%.2f", s)) J/mol/K")
        println("  G  = $(@sprintf("%.2f", g)) kJ/mol")
    catch e
        println("  Hierarchical selection failed: $(e)")
    end
    
    # Process all species
    println("\nProcessing all species with 4 properties per log plot...")
    
    # Get all species
    all_species = get_all_species(db)
    println("Found $(length(all_species)) species to process")
    
    # Process each species
    success_count = 0
    success_species = String[]
    
    for (i, species) in enumerate(all_species)
        print("Processing $(i)/$(length(all_species)): $(species)...")
        success = plot_all_properties(db, species)
        if success
            success_count += 1
            push!(success_species, species)
            println(" Done!")
        else
            println(" Failed!")
        end
    end
    
    # Create summary report
    summary_df = create_summary_report(db, all_species, success_count, success_species)
    
    println("\nSuccessfully processed $(success_count) out of $(length(all_species)) species")
    println("Plots saved to the 'plots' directory.")
    println("Tabular data saved to 'output/tables' directory.")
    println("Source documentation saved to 'output/docs' directory.")
    
    println("\nThis implementation reads data directly from original sources")
    println("using the hierarchical approach, prioritizing data from the")
    println("most respected sources, with log temperature scales and uncertainty shading.")
    
    # Print summary of sources used
    source_counts = Dict()
    for src in summary_df.Source
        source_counts[src] = get(source_counts, src, 0) + 1
    end
    
    println("\nSource usage summary:")
    for (src, count) in sort(collect(source_counts), by=x->x[2], rev=true)
        println("  $(rpad(src, 15)) : $(count) species")
    end
end

# Run the main function
main()