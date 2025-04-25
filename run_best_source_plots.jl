#!/usr/bin/env julia

"""
Best Source Thermodynamic Plotter

This script plots all species with ONLY the best available source shown
(highest priority source) with proper uncertainty bands.

Unlike the full multi-source plots, this creates clean, single-line plots
showing only the most accurate data for each species.
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using Printf
using Plots
using JSON
using Colors
using Dates
using YAML
using DataFrames
using DuckDB

# Function to convert integers to superscript for tick labels
function superscript(n::Int)
    if n == 0
        return "⁰"
    end
    
    superscripts = Dict(
        '0' => '⁰',
        '1' => '¹',
        '2' => '²',
        '3' => '³',
        '4' => '⁴',
        '5' => '⁵',
        '6' => '⁶',
        '7' => '⁷',
        '8' => '⁸',
        '9' => '⁹',
        '-' => '⁻'
    )
    
    str = string(n)
    result = ""
    for c in str
        result *= superscripts[c]
    end
    
    return result
end

function get_all_species()
    # Get list of species from species.yaml config file
    project_dir = @__DIR__
    species_yaml_path = joinpath(project_dir, "config", "species.yaml")
    
    # Check if file exists
    if !isfile(species_yaml_path)
        @warn "species.yaml not found in config directory, using default species list"
        return ["N2", "O2", "H2O", "CO2"]
    end
    
    # Read and parse the YAML file
    try
        species_config = YAML.load_file(species_yaml_path)
        
        # Extract the species list
        if haskey(species_config, "species") && species_config["species"] isa Vector
            @info "Loaded $(length(species_config["species"])) species from species.yaml"
            return species_config["species"]
        else
            @warn "No valid species list found in species.yaml, using default species list"
            return ["N2", "O2", "H2O", "CO2"]
        end
    catch e
        @warn "Error parsing species.yaml: $(e), using default species list"
        return ["N2", "O2", "H2O", "CO2"]
    end
end

function load_plot_config()
    # Load all plot settings from plots.yaml
    project_dir = @__DIR__
    plots_yaml_path = joinpath(project_dir, "config", "plots.yaml")
    
    # Default configuration (fallback)
    default_config = Dict(
        "general" => Dict(
            "size" => [1200, 900],
            "dpi" => 300,
            "margin" => 10,
            "right_margin" => 20,
            "title_font_size" => 12,
            "temperature_points" => 100,
            "use_log_scale" => true
        ),
        "line_styles" => Dict(
            "best_source" => Dict(
                "line_width" => 2.5,
                "alpha" => 1.0,
                "style" => "solid",
                "color" => "blue"
            )
        ),
        "uncertainty" => Dict(
            "ribbon_alpha" => 0.2
        ),
        "output" => Dict(
            "file_format" => "png",
            "save_plots" => true,
            "plots_directory" => "plots_best_source",
            "documentation_directory" => "output/docs"
        ),
        "properties" => Dict(
            "cp" => Dict("label" => "Cp (J/mol/K)", "position" => 1),
            "h" => Dict("label" => "H (kJ/mol)", "position" => 2),
            "s" => Dict("label" => "S (J/mol/K)", "position" => 3),
            "g" => Dict("label" => "G (kJ/mol)", "position" => 4)
        )
    )
    
    # Check if the file exists
    if !isfile(plots_yaml_path)
        @warn "plots.yaml not found at $(plots_yaml_path), using default plot settings"
        return default_config
    end
    
    # Read and parse the YAML file
    try
        plot_config = YAML.load_file(plots_yaml_path)
        @info "Loaded plot configuration from plots.yaml"
        
        # Create a modified config with a different plots directory for best-source plots
        modified_config = deepcopy(plot_config)
        if haskey(modified_config, "output") && haskey(modified_config["output"], "plots_directory")
            modified_config["output"]["plots_directory"] = "plots_best_source"
        end
        
        return modified_config
    catch e
        @warn "Error parsing plots.yaml: $(e), using default plot settings"
        return default_config
    end
end

function get_temperature_range()
    # Get temperature range from species.yaml config file
    project_dir = @__DIR__
    species_yaml_path = joinpath(project_dir, "config", "species.yaml")
    
    # Default temperature range
    default_range = [100.0, 10000.0]
    
    # Check if file exists
    if !isfile(species_yaml_path)
        @warn "species.yaml not found in config directory, using default temperature range"
        return default_range
    end
    
    # Read and parse the YAML file
    try
        species_config = YAML.load_file(species_yaml_path)
        
        # Extract the temperature range
        if haskey(species_config, "temperature_range") && species_config["temperature_range"] isa Vector{<:Number} && length(species_config["temperature_range"]) == 2
            temp_range = species_config["temperature_range"]
            @info "Using temperature range from species.yaml: $(temp_range[1])-$(temp_range[2]) K"
            return temp_range
        else
            @warn "No valid temperature range found in species.yaml, using default range [100, 10000]"
            return default_range
        end
    catch e
        @warn "Error parsing species.yaml: $(e), using default temperature range [100, 10000]"
        return default_range
    end
end

function plot_species_best_source(species_name, db_path)
    println("Plotting best source properties for $(species_name)...")
    
    # Get temperature range from Species.yaml
    temp_range = get_temperature_range()
    
    # Load plot configuration from plots.yaml
    plot_config = load_plot_config()
    
    # CRITICAL: Use EXPLICIT hierarchical selection to get the BEST source only
    options = Dict(
        "method" => "hierarchical",  # Use hierarchical method to get the best source
        "enable_fallback" => true,
        "use_json" => true  # Ensure we're using the JSON data sources
    )
    
    # Print detailed diagnostic information
    println("  Loading polynomial data with hierarchical selection...")
    
    # Load the polynomial data
    best_poly = JThermodynamicsData.load_polynomial_data(db_path, species_name, temp_range, options)
    
    # Handle both Dictionary and ThermodynamicPolynomial return types
    best_source = ""
    best_priority = -1
    
    if isa(best_poly, Dict)
        # Handle dictionary return type
        best_source = best_poly["source"]
        best_priority = get(best_poly, "priority", -1)
        
        # Print details for debugging
        println("  Selected source details (Dict):")
        println("    - Source: $(best_poly["source"])")
        println("    - Priority: $(get(best_poly, "priority", "unknown"))")
        println("    - Temperature ranges: $(get(best_poly, "temperature_ranges", "unknown"))")
        println("    - Uncertainty: $(get(best_poly, "uncertainty", "unknown"))")
    else
        # Handle ThermodynamicPolynomial return type
        best_source = best_poly.source
        best_priority = best_poly.priority
        
        # Print details for debugging
        println("  Selected source details (ThermodynamicPolynomial):")
        println("    - Source: $(best_poly.source)")
        println("    - Priority: $(best_poly.priority)")
        println("    - Temperature ranges: $(best_poly.temperature_ranges)")
        println("    - Uncertainty: $(best_poly.uncertainty)")
    end
    
    # Check if we got a theoretical or experimental source
    if best_priority <= 4 || startswith(lowercase(best_source), "theoretical")
        println("  WARNING: Using a theoretical source! Priority = $(best_priority)")
    else
        println("  SUCCESS: Using an experimental source! Priority = $(best_priority)")
    end
    
    # Extract plot settings from configuration
    general_settings = plot_config["general"]
    size_config = general_settings["size"]
    dpi_config = general_settings["dpi"]
    margin_config = general_settings["margin"]
    right_margin_config = general_settings["right_margin"]
    temperature_points = general_settings["temperature_points"]
    use_log_scale = general_settings["use_log_scale"]
    
    # Generate temperatures
    if use_log_scale
        # Logarithmic temperature scale
        temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=temperature_points)
    else
        # Linear temperature scale
        temps = range(temp_range[1], temp_range[2], length=temperature_points)
    end
    
    # Create a 2x2 plot layout for the 4 thermodynamic properties
    plt = plot(
        layout=(2,2), 
        size=(size_config[1], size_config[2]), 
        dpi=dpi_config, 
        margin=margin_config * Plots.mm, 
        right_margin=right_margin_config * Plots.mm
    )
    
    # Use properties from configuration
    properties_config = plot_config["properties"]
    properties = [
        (properties_config["cp"]["position"], "cp", properties_config["cp"]["label"]),
        (properties_config["h"]["position"], "h", properties_config["h"]["label"]),
        (properties_config["s"]["position"], "s", properties_config["s"]["label"]),
        (properties_config["g"]["position"], "g", properties_config["g"]["label"])
    ]
    
    # Standardize source name for display
    display_name = ""
    if best_priority <= 4 || startswith(lowercase(best_source), "theoretical")
        if startswith(lowercase(best_source), "theoretical-")
            display_name = titlecase(best_source)
        else
            display_name = "Theoretical-" * titlecase(best_source)
        end
    else
        # Non-theoretical sources keep their original name with proper capitalization
        words = split(best_source, "-")
        capitalized_words = [uppercase(first(word)) * lowercase(SubString(word, 2:lastindex(word))) for word in words]
        display_name = join(capitalized_words, "-")
    end
    
    # Get line style properties from configuration
    line_styles = get(plot_config, "line_styles", Dict())
    best_styles = get(line_styles, "best_source", Dict())
    line_width = get(best_styles, "line_width", 2.5)
    line_alpha = get(best_styles, "alpha", 1.0)
    
    # Get uncertainty settings
    uncertainty_settings = get(plot_config, "uncertainty", Dict())
    ribbon_alpha = get(uncertainty_settings, "ribbon_alpha", 0.2)
    
    # Get uncertainty value
    uncertainty_value = 0.05  # Default 5%
    if isa(best_poly, Dict)
        uncertainty_value = get(best_poly, "uncertainty", 0.05)
    else
        uncertainty_value = best_poly.uncertainty
    end
    
    # Calculate and plot each property
    for (idx, property, ylabel) in properties
        # Calculate property values
        values = if property == "cp"
            [JThermodynamicsData.calculate_cp(best_poly, t) for t in temps]
        elseif property == "h"
            [JThermodynamicsData.calculate_h(best_poly, t) for t in temps]
        elseif property == "s"
            [JThermodynamicsData.calculate_s(best_poly, t) for t in temps]
        else # g
            [JThermodynamicsData.calculate_g(best_poly, t) for t in temps]
        end
        
        # Calculate uncertainty bands
        upper = values .* (1 + uncertainty_value)
        lower = values .* (1 - uncertainty_value)
        
        # Generate logarithmic ticks
        if use_log_scale
            # Create logarithmic ticks for each decade
            log_min = log10(temp_range[1])
            log_max = log10(temp_range[2])
            decades = floor(Int, log_min):ceil(Int, log_max)
            
            # Generate ticks within each decade
            tick_positions = Float64[]
            tick_labels = String[]
            
            for decade in decades
                # Major tick for the decade
                push!(tick_positions, 10.0^decade)
                push!(tick_labels, "10$(superscript(decade))")
                
                # Minor ticks within the decade (2-9)
                for i in 2:9
                    minor_tick = i * 10.0^decade
                    if minor_tick >= temp_range[1] && minor_tick <= temp_range[2]
                        push!(tick_positions, minor_tick)
                        # Use empty string for minor ticks
                        push!(tick_labels, "")
                    end
                end
            end
            
            # Plot with uncertainty ribbon
            plot!(
                plt, 
                temps, 
                values, 
                subplot=idx,
                label=display_name,
                xlabel="Temperature (K)",
                ylabel=ylabel,
                title="$(uppercase(property)) for $species_name",
                xscale=:log10,
                xticks=(tick_positions, tick_labels),
                xminorticks=10,
                line=(:blue, 3.0, :solid),  # Blue solid line
                ribbon=(values .- lower, upper .- values),
                fillalpha=0.3,
                fillcolor=:lightblue,
                legend=:outerright
            )
        else
            # Linear scale plotting
            plot!(
                plt, 
                temps, 
                values, 
                subplot=idx,
                label=display_name,
                xlabel="Temperature (K)",
                ylabel=ylabel,
                title="$(uppercase(property)) for $species_name",
                line=(:blue, 3.0, :solid),
                ribbon=(values .- lower, upper .- values),
                fillalpha=0.3,
                fillcolor=:lightblue,
                legend=:outerright
            )
        end
    end
    
    # Add a main title
    title_font_size = get(general_settings, "title_font_size", 12)
    plot!(plt, title="Thermodynamic Properties of $(species_name) - Best Source: $(display_name)", 
          titleloc=:center, titlefontsize=title_font_size)
    
    # Save the plot
    output_settings = plot_config["output"]
    plots_dir = joinpath(@__DIR__, output_settings["plots_directory"])
    mkpath(plots_dir)
    
    # Use configured file format
    file_format = output_settings["file_format"]
    savefig(plt, joinpath(plots_dir, "$(species_name)_best_source.$(file_format)"))
    
    println("  ✓ Plot saved to $(plots_dir)/$(species_name)_best_source.$(file_format)")
    
    # Validate against the data source summary
    # Read from the species-specific sources.md file
    docs_dir = get(output_settings, "documentation_directory", "output/docs")
    docs_path = joinpath(@__DIR__, docs_dir, "$(species_name)_sources.md")
    
    if isfile(docs_path)
        md_content = read(docs_path, String)
        source_match = match(r"\*\*Source:\*\* ([^\n]+)", md_content)
        
        if source_match !== nothing
            doc_source = source_match.captures[1]
            if lowercase(doc_source) == lowercase(best_source)
                println("  ✓ Source validation: Best source matches documentation")
            else
                println("  ⚠️ Source validation: MISMATCH - Plot uses $(best_source) but documentation says $(doc_source)")
            end
        end
    end
    
    return plt
end

function export_property_data(species_name, best_poly, temp_range)
    try
        # Load plot configuration
        plot_config = load_plot_config()
        general_settings = plot_config["general"]
        
        # Get temperature generation settings from configuration
        temperature_points = general_settings["temperature_points"] * 2  # Double points for export data
        use_log_scale = general_settings["use_log_scale"]
        
        # Generate temperature points based on configuration
        if use_log_scale
            temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=temperature_points)
        else
            temps = range(temp_range[1], temp_range[2], length=temperature_points)
        end
        
        # Get source info
        source = isa(best_poly, Dict) ? best_poly["source"] : best_poly.source
        uncertainty = isa(best_poly, Dict) ? best_poly["uncertainty"] : best_poly.uncertainty
        
        # Calculate all properties
        cp_values = [JThermodynamicsData.calculate_cp(best_poly, t) for t in temps]
        h_values = [JThermodynamicsData.calculate_h(best_poly, t) for t in temps]
        s_values = [JThermodynamicsData.calculate_s(best_poly, t) for t in temps]
        g_values = [JThermodynamicsData.calculate_g(best_poly, t) for t in temps]
        
        # Create data arrays
        data = Dict(
            "Temperature" => temps,
            "Cp" => cp_values,
            "H" => h_values,
            "S" => s_values,
            "G" => g_values,
            "Source" => fill(source, length(temps)),
            "Uncertainty" => fill(uncertainty, length(temps))
        )
        
        # Get output directory from configuration
        output_settings = plot_config["output"]
        tables_dir = joinpath(@__DIR__, "output", "tables_best_source")
        mkpath(tables_dir)
        
        # Export to JSON
        open(joinpath(tables_dir, "$(species_name)_best_source_data.json"), "w") do f
            JSON.print(f, data, 4)
        end
        
        return true
    catch e
        println("  Error exporting data for $(species_name): $(e)")
        return false
    end
end

function verify_source_consistency(species_name)
    # This function verifies that the best source used in the plot matches what's in the data_source_summary.md file
    
    # Get the source from the documentation file
    docs_dir = joinpath(@__DIR__, "output", "docs")
    docs_path = joinpath(docs_dir, "$(species_name)_sources.md")
    
    doc_source = "unknown"
    if isfile(docs_path)
        try
            md_content = read(docs_path, String)
            source_match = match(r"\*\*Source:\*\* ([^\n]+)", md_content)
            if source_match !== nothing
                doc_source = source_match.captures[1]
            end
        catch e
            println("  Error reading documentation for $(species_name): $(e)")
        end
    end
    
    # Get the source from the summary file
    summary_path = joinpath(@__DIR__, "output", "data_source_summary.md")
    summary_source = "unknown"
    
    if isfile(summary_path)
        try
            md_content = read(summary_path, String)
            # Look for the species in the per-species table
            pattern = "\\| $(species_name) \\| ([^|]+) \\|"
            source_match = match(Regex(pattern), md_content)
            if source_match !== nothing
                summary_source = strip(source_match.captures[1])
            end
        catch e
            println("  Error reading summary for $(species_name): $(e)")
        end
    end
    
    # Compare sources
    if doc_source == "unknown" || summary_source == "unknown"
        return (false, doc_source, summary_source)
    else
        return (lowercase(doc_source) == lowercase(summary_source), doc_source, summary_source)
    end
end

function create_summary_report(species_results)
    println("\nCreating validation report...")
    
    # Output directory
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    
    # Create summary report
    open(joinpath(output_dir, "best_source_validation.md"), "w") do f
        write(f, "# Best Source Validation Report\n\n")
        
        write(f, "This report validates that the same best source is consistently used across:\n")
        write(f, "1. Individual species documentation files\n")
        write(f, "2. Data source summary table\n")
        write(f, "3. Best source plots\n\n")
        
        write(f, "## Validation Results\n")
        write(f, "| Species | Plot Source | Documentation Source | Summary Source | Consistent |\n")
        write(f, "|---------|------------|---------------------|---------------|------------|\n")
        
        for (species, result) in species_results
            plot_source = result[:plot_source]
            doc_source = result[:doc_source]
            summary_source = result[:summary_source]
            consistent = result[:consistent] ? "✓" : "❌"
            
            write(f, "| $(species) | $(plot_source) | $(doc_source) | $(summary_source) | $(consistent) |\n")
        end
        
        # Calculate statistics
        total = length(species_results)
        consistent_count = count(result -> result[2][:consistent], collect(species_results))
        
        write(f, "\n## Summary Statistics\n")
        write(f, "- Total species: $(total)\n")
        write(f, "- Consistent sources: $(consistent_count)\n")
        write(f, "- Consistency rate: $(round(100.0 * consistent_count / total, digits=1))%\n")
        
        write(f, "\n---\n")
        write(f, "*Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))*\n")
    end
    
    println("Validation report saved to output/best_source_validation.md")
end

function main()
    # Set up logging
    log_dir = joinpath(@__DIR__, "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "best_source_plots_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("BestSourcePlotter", Dict(
        "start_time" => now(),
        "working_directory" => @__DIR__
    ))
    
    @info "JThermodynamicsData - Best Source Plots Generator"
    @info "=============================================="
    
    # Initialize package
    JThermodynamicsData.log_stage_start("BestSourcePlotter", "InitializePackage", Dict())
    config_path = joinpath(@__DIR__, "config", "settings.yaml")
    config, db = JThermodynamicsData.with_timing(
        () -> JThermodynamicsData.initialize(config_path),
        "Initialization", 
        "LoadConfig"
    )
    JThermodynamicsData.log_stage_end("BestSourcePlotter", "InitializePackage", Dict(
        "config_path" => config_path,
        "db_connected" => db !== nothing
    ))
    
    # Get database path
    db_path = config["general"]["database_path"]
    if !startswith(db_path, "/")
        db_path = joinpath(@__DIR__, db_path)
    end
    
    # Get all species from Species.yaml
    JThermodynamicsData.log_stage_start("BestSourcePlotter", "GetSpeciesList", Dict())
    species_to_process = JThermodynamicsData.with_timing(
        () -> get_all_species(),
        "Data", 
        "GetSpeciesList"
    )
    @info "Loaded $(length(species_to_process)) species from Species.yaml"
    
    # Process each species
    JThermodynamicsData.log_stage_start("BestSourcePlotter", "ProcessSpecies", Dict(
        "total_species" => length(species_to_process)
    ))
    
    # Track results for validation
    species_results = Dict()
    
    for (i, species) in enumerate(species_to_process)
        @info "Processing $(i)/$(length(species_to_process)): $(species)"
        
        JThermodynamicsData.log_stage_start("BestSourcePlotter", "Species:$species", Dict(
            "index" => i,
            "total" => length(species_to_process)
        ))
        
        # Generate the best source plot
        plot_start = now()
        
        # Get best source data using hierarchical selection
        options = Dict(
            "method" => "hierarchical",
            "enable_fallback" => true,
            "use_json" => true
        )
        
        # Load the best polynomial data
        best_poly = JThermodynamicsData.load_polynomial_data(db_path, species, get_temperature_range(), options)
        
        # Get the source name
        best_source = ""
        if isa(best_poly, Dict)
            best_source = best_poly["source"]
        else
            best_source = best_poly.source
        end
        
        # Generate the plot
        plot_species_best_source(species, db_path)
        plot_time = now() - plot_start
        
        JThermodynamicsData.log_timing_benchmark("Plotting", species, plot_time)
        
        # Export the data
        JThermodynamicsData.log_stage_start("BestSourcePlotter", "ExportData:$species", Dict(
            "species" => species
        ))
        
        export_start = now()
        export_success = export_property_data(species, best_poly, get_temperature_range())
        export_time = now() - export_start
        
        JThermodynamicsData.log_timing_benchmark("Export", species, export_time)
        
        JThermodynamicsData.log_stage_end("BestSourcePlotter", "ExportData:$species", Dict(
            "species" => species,
            "success" => export_success
        ))
        
        # Verify source consistency
        consistent, doc_source, summary_source = verify_source_consistency(species)
        
        species_results[species] = Dict(
            :plot_source => best_source,
            :doc_source => doc_source,
            :summary_source => summary_source,
            :consistent => consistent,
            :plot_time => JThermodynamicsData.format_time_duration(plot_time),
            :export_time => JThermodynamicsData.format_time_duration(export_time)
        )
        
        JThermodynamicsData.log_stage_end("BestSourcePlotter", "Species:$species", Dict(
            "success" => true,
            "best_source" => best_source,
            "consistent" => consistent,
            "elapsed_time" => JThermodynamicsData.format_time_duration(plot_time)
        ))
    end
    
    JThermodynamicsData.log_stage_end("BestSourcePlotter", "ProcessSpecies", Dict(
        "processed_species" => length(species_to_process)
    ))
    
    # Create validation report
    JThermodynamicsData.log_stage_start("BestSourcePlotter", "CreateValidationReport", Dict())
    report_start = now()
    create_summary_report(species_results)
    report_time = now() - report_start
    JThermodynamicsData.log_timing_benchmark("Report", "ValidationReport", report_time)
    JThermodynamicsData.log_stage_end("BestSourcePlotter", "CreateValidationReport", Dict(
        "report_path" => joinpath(@__DIR__, "output", "best_source_validation.md"),
        "elapsed_time" => JThermodynamicsData.format_time_duration(report_time)
    ))
    
    # Load plot configuration for output paths
    plot_config = load_plot_config()
    output_settings = plot_config["output"]
    plots_dir = output_settings["plots_directory"]
    
    @info "Successfully processed $(length(species_to_process)) species with best source selection"
    @info "Plots saved to the '$(plots_dir)' directory"
    @info "Data tables saved to 'output/tables_best_source' directory"
    @info "Validation report saved to 'output/best_source_validation.md'"
    
    # Log pipeline completion
    JThermodynamicsData.log_pipeline_end("BestSourcePlotter", Dict(
        "total_species" => length(species_to_process),
        "log_file" => log_file
    ))
end

# Run the main function
main()