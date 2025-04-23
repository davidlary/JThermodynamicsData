#!/usr/bin/env julia

# All Species Multi-Source Thermodynamic Plotter
# This script plots all species with all available sources shown
# using hierarchical selection for the best source

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using Printf
using Plots
using JSON
using Colors
using Dates

function get_all_species()
    # Get list of all species from JSON storage
    return JThermodynamicsData.list_available_species()
end

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
    db_path = joinpath(@__DIR__, "data", "thermodynamics.duckdb")
    
    # Load best source data (will be displayed prominently)
    best_poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, Dict("method" => "hierarchical"))
    best_source = best_poly.source
    
    println("  Selected source: $(best_source)")
    
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
    
    # Standardize theoretical sources group for consistent legend formatting
    theoretical_sources = [
        "theoretical", 
        "THEORETICAL_GROUP_CONTRIBUTION", 
        "THEORETICAL_STATISTICAL_THERMO",
        "THEORETICAL_MACHINE_LEARNING", 
        "THEORETICAL_QUANTUM_CHEMISTRY",
        "THEORETICAL_ENSEMBLE",
        "THEORETICAL_FALLBACK"
    ]
    
    # Map sources to their standardized display name
    source_display_names = Dict()
    for (src_name, _, _) in sources
        # Standardize theoretical source names for consistent legends
        if startswith(lowercase(src_name), "theoretical")
            if lowercase(src_name) == "theoretical"
                source_display_names[src_name] = "theoretical"
            elseif lowercase(src_name) == "theoretical_group_contribution"
                source_display_names[src_name] = "group-contribution"
            elseif lowercase(src_name) == "theoretical_statistical_thermo"
                source_display_names[src_name] = "stat-thermo"
            elseif lowercase(src_name) == "theoretical_machine_learning"
                source_display_names[src_name] = "machine-learning"
            elseif lowercase(src_name) == "theoretical_quantum_chemistry"
                source_display_names[src_name] = "quantum-statistical"
            elseif lowercase(src_name) == "theoretical_ensemble"
                source_display_names[src_name] = "benson-group"
            else
                source_display_names[src_name] = "theoretical"
            end
        else
            # Non-theoretical sources keep their original name
            source_display_names[src_name] = src_name
        end
    end
    
    # Ensure all species use consistent color scheme
    # Blue for best source, gray scale for others, consistent colors for theoretical sources
    source_colors = Dict()
    line_styles = Dict()
    
    # Define a consistent color palette for theoretical sources
    theoretical_colors = Dict(
        "theoretical" => RGB(0.7, 0.7, 0.7),
        "group-contribution" => RGB(0.6, 0.6, 0.6),
        "stat-thermo" => RGB(0.5, 0.5, 0.5),
        "machine-learning" => RGB(0.45, 0.45, 0.45),
        "quantum-statistical" => RGB(0.4, 0.4, 0.4),
        "benson-group" => RGB(0.35, 0.35, 0.35)
    )
    
    # Assign colors and line styles
    for (i, (src_name, priority, _)) in enumerate(sources)
        display_name = source_display_names[src_name]
        
        if src_name == best_source
            source_colors[src_name] = :blue
            line_styles[src_name] = :solid
        elseif haskey(theoretical_colors, display_name)
            # Use predefined theoretical colors
            source_colors[src_name] = theoretical_colors[display_name]
            line_styles[src_name] = :dash
        else
            # Non-theoretical sources get color by priority
            priority_level = max(1, min(priority, 13))
            color_value = 0.8 - (priority_level / 20.0)  # Higher priority = darker color
            source_colors[src_name] = RGB(color_value, color_value, color_value)
            line_styles[src_name] = :dash
        end
    end
    
    # Plot each property for each source
    for (source_name, priority, reliability) in sources
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
                    plot!(
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
    plot!(plt, title="Thermodynamic Properties for $(species_name)", 
          titleloc=:center, titlefontsize=12)
    
    # Save the plot
    plots_dir = joinpath(@__DIR__, "plots")
    mkpath(plots_dir)
    savefig(plt, joinpath(plots_dir, "$(species_name)_all_sources.png"))
    
    # Create source documentation in Markdown
    docs_dir = joinpath(@__DIR__, "output", "docs")
    mkpath(docs_dir)
    
    open(joinpath(docs_dir, "$(species_name)_sources.md"), "w") do f
        write(f, "# Thermodynamic Data Source for $(species_name)\n\n")
        write(f, "## Selected Source\n")
        write(f, "- **Source:** $(best_poly.source)\n")
        write(f, "- **Priority:** $(best_poly.priority)\n")
        write(f, "- **Reliability score:** $(best_poly.reliability_score)\n")
        write(f, "- **Uncertainty:** $(best_poly.uncertainty * 100)%\n")
        write(f, "- **Temperature range:** $(temp_range[1]) - $(temp_range[2]) K\n")
        write(f, "- **Polynomial type:** $(best_poly.polynomial_type)\n\n")
        
        write(f, "## Coefficients\n")
        write(f, "### Low temperature range ($(best_poly.temperature_ranges[1][1]) - $(best_poly.temperature_ranges[1][2]) K)\n")
        write(f, "```\n$(join([@sprintf("%.8e", c) for c in best_poly.coefficients[1]], ", "))\n```\n\n")
        write(f, "### High temperature range ($(best_poly.temperature_ranges[2][1]) - $(best_poly.temperature_ranges[2][2]) K)\n")
        write(f, "```\n$(join([@sprintf("%.8e", c) for c in best_poly.coefficients[2]], ", "))\n```\n")
        
        # Add information about all available sources
        write(f, "\n## All Available Sources\n")
        write(f, "| Source | Priority | Reliability |\n")
        write(f, "|--------|----------|-------------|\n")
        for (src_name, priority, reliability) in sources
            # Use standardized display name for documentation
            display_name = source_display_names[src_name]
            write(f, "| $(display_name) | $(priority) | $(reliability) |\n")
        end
    end
    
    println("  Plot saved to plots/$(species_name)_all_sources.png")
    println("  Documentation saved to output/docs/$(species_name)_sources.md")
    return true
end

function export_property_data(species_name, temp_range=[200.0, 6000.0])
    try
        # Get database path
        db_path = joinpath(@__DIR__, "data", "thermodynamics.duckdb")
        
        # Load data using hierarchical selection
        poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, Dict("method" => "hierarchical"))
        
        # Generate temperature points (using logarithmic scale)
        temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=200)
        
        # Calculate all properties
        cp_values = [JThermodynamicsData.calculate_cp(poly, t) for t in temps]
        h_values = [JThermodynamicsData.calculate_h(poly, t) for t in temps]
        s_values = [JThermodynamicsData.calculate_s(poly, t) for t in temps]
        g_values = [JThermodynamicsData.calculate_g(poly, t) for t in temps]
        
        # Create data arrays
        data = Dict(
            "Temperature" => temps,
            "Cp" => cp_values,
            "H" => h_values,
            "S" => s_values,
            "G" => g_values,
            "Source" => fill(poly.source, length(temps)),
            "Uncertainty" => fill(poly.uncertainty, length(temps))
        )
        
        # Create output directory
        tables_dir = joinpath(@__DIR__, "output", "tables")
        mkpath(tables_dir)
        
        # Export to JSON
        open(joinpath(tables_dir, "$(species_name)_data.json"), "w") do f
            JSON.print(f, data, 4)
        end
        
        # Also export the polynomial coefficients
        coeff_data = Dict(
            "species" => species_name,
            "source" => poly.source,
            "priority" => poly.priority,
            "reliability_score" => poly.reliability_score,
            "polynomial_type" => poly.polynomial_type,
            "temperature_ranges" => poly.temperature_ranges,
            "coefficients" => poly.coefficients,
            "uncertainty" => poly.uncertainty
        )
        
        open(joinpath(tables_dir, "$(species_name)_coefficients.json"), "w") do f
            JSON.print(f, coeff_data, 4)
        end
        
        return true
    catch e
        println("  Error exporting data for $(species_name): $(e)")
        return false
    end
end

function create_summary_report(all_species, success_species)
    println("\nCreating summary report...")
    
    # Output directory
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    
    # Create source statistics
    source_counts = Dict()
    total_sources = 0
    
    for species in success_species
        sources = JThermodynamicsData.list_available_sources(species)
        for (source_name, _, _) in sources
            source_counts[source_name] = get(source_counts, source_name, 0) + 1
            total_sources += 1
        end
    end
    
    # Create main summary report
    open(joinpath(output_dir, "summary_report.md"), "w") do f
        write(f, "# Thermodynamic Data Summary Report\n\n")
        
        write(f, "## Overview\n")
        write(f, "- Total species: $(length(all_species))\n")
        write(f, "- Successfully processed: $(length(success_species))\n")
        write(f, "- Total data sources used: $(total_sources)\n\n")
        
        write(f, "## Source Statistics\n")
        write(f, "| Source | Count | Percentage |\n")
        write(f, "|--------|-------|------------|\n")
        
        for (source, count) in sort(collect(source_counts), by=x->x[2], rev=true)
            percentage = round(count / total_sources * 100, digits=1)
            write(f, "| $(source) | $(count) | $(percentage)% |\n")
        end
        
        write(f, "\n## Species with Multiple Sources\n")
        write(f, "| Species | Sources | Best Source |\n")
        write(f, "|---------|---------|-------------|\n")
        
        # Find species with multiple sources
        for species in success_species
            sources = JThermodynamicsData.list_available_sources(species)
            if length(sources) > 1
                # Get best source
                db_path = joinpath(@__DIR__, "data", "thermodynamics.duckdb")
                try
                    poly = JThermodynamicsData.load_polynomial_data(db_path, species, [200.0, 6000.0], Dict("method" => "hierarchical"))
                    best_source = poly.source
                    write(f, "| $(species) | $(length(sources)) | $(best_source) |\n")
                catch e
                    write(f, "| $(species) | $(length(sources)) | error |\n")
                end
            end
        end
        
        write(f, "\n## All Species List\n")
        write(f, "| Species | Processed | Sources |\n")
        write(f, "|---------|-----------|--------|\n")
        
        for species in sort(all_species)
            processed = species in success_species ? "✓" : "❌"
            source_count = 0
            if species in success_species
                source_count = length(JThermodynamicsData.list_available_sources(species))
            end
            write(f, "| $(species) | $(processed) | $(source_count) |\n")
        end
    end
    
    println("Summary report saved to output/summary_report.md")
end

function main()
    # Set up logging
    log_dir = joinpath(@__DIR__, "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "all_species_plots_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("AllSpeciesPlotter", Dict(
        "start_time" => now(),
        "working_directory" => @__DIR__
    ))
    
    @info "JThermodynamicsData - All Species Plots Generator"
    @info "=============================================="
    
    # Initialize package
    JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "InitializePackage", Dict())
    config_path = joinpath(@__DIR__, "config", "settings.yaml")
    config, db = JThermodynamicsData.with_timing(
        () -> JThermodynamicsData.initialize(config_path),
        "Initialization", 
        "LoadConfig"
    )
    JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "InitializePackage", Dict(
        "config_path" => config_path,
        "db_connected" => db !== nothing
    ))
    
    # Get all species
    JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "GetSpeciesList", Dict())
    all_species = JThermodynamicsData.with_timing(
        () -> get_all_species(),
        "Data", 
        "GetSpeciesList"
    )
    @info "Found $(length(all_species)) species"
    JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "GetSpeciesList", Dict(
        "species_count" => length(all_species)
    ))
    
    # Process each species
    JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "ProcessSpecies", Dict(
        "total_species" => length(all_species)
    ))
    
    success_species = String[]
    species_results = Dict()
    
    for (i, species) in enumerate(all_species)
        @info "Processing $(i)/$(length(all_species)): $(species)"
        
        JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "Species:$species", Dict(
            "index" => i,
            "total" => length(all_species)
        ))
        
        # Plot the species with all sources
        plot_start = now()
        success = plot_species_with_all_sources(species)
        plot_time = now() - plot_start
        
        JThermodynamicsData.log_timing_benchmark("Plotting", species, plot_time)
        
        if success
            push!(success_species, species)
            
            # Also export data tables
            JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "ExportData:$species", Dict(
                "species" => species
            ))
            
            export_start = now()
            export_success = export_property_data(species)
            export_time = now() - export_start
            
            JThermodynamicsData.log_timing_benchmark("Export", species, export_time)
            
            JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "ExportData:$species", Dict(
                "species" => species,
                "success" => export_success
            ))
            
            species_results[species] = Dict(
                "plot_success" => true,
                "export_success" => export_success,
                "plot_time" => JThermodynamicsData.format_time_duration(plot_time),
                "export_time" => JThermodynamicsData.format_time_duration(export_time)
            )
        else
            species_results[species] = Dict(
                "plot_success" => false
            )
        end
        
        JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "Species:$species", Dict(
            "success" => success,
            "elapsed_time" => JThermodynamicsData.format_time_duration(plot_time)
        ))
    end
    
    JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "ProcessSpecies", Dict(
        "processed_species" => length(all_species),
        "successful_species" => length(success_species),
        "success_rate" => "$(round(length(success_species) / length(all_species) * 100, digits=1))%"
    ))
    
    # Create summary report
    JThermodynamicsData.log_stage_start("AllSpeciesPlotter", "CreateSummaryReport", Dict())
    summary_start = now()
    create_summary_report(all_species, success_species)
    summary_time = now() - summary_start
    JThermodynamicsData.log_timing_benchmark("Report", "SummaryReport", summary_time)
    JThermodynamicsData.log_stage_end("AllSpeciesPlotter", "CreateSummaryReport", Dict(
        "report_path" => joinpath(@__DIR__, "output", "summary_report.md"),
        "elapsed_time" => JThermodynamicsData.format_time_duration(summary_time)
    ))
    
    @info "Successfully processed $(length(success_species)) out of $(length(all_species)) species"
    @info "Plots saved to the 'plots' directory"
    @info "Data tables saved to 'output/tables' directory"
    @info "Source documentation saved to 'output/docs' directory"
    @info "Main summary report saved to 'output/summary_report.md'"
    @info "This implementation uses the hierarchical selection system, showing all available sources with the best one highlighted."
    
    # Log pipeline completion
    JThermodynamicsData.log_pipeline_end("AllSpeciesPlotter", Dict(
        "total_species" => length(all_species),
        "successful_species" => length(success_species),
        "success_rate" => "$(round(length(success_species) / length(all_species) * 100, digits=1))%",
        "log_file" => log_file
    ))
end

# Run the main function
main()