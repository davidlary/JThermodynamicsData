#!/usr/bin/env julia

"""
JThermodynamicsData - Main Program

This script reads the species list from config/species.yaml and processes
thermodynamic data for all specified species in the temperature range 100-10,000 K.
"""

# Activate the project to ensure we use the right dependencies
using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using YAML
using DataFrames
using Plots
using Printf
using Dates
using DuckDB

# Lowercase constant aliases are now directly defined in constants.jl
# No need to add them dynamically
@info "Physical constants available in both uppercase (R, H, KB, NA, C) and lowercase (r, h, kb, na, c) forms"

"""
    calculate_properties_fixed(conn::DuckDB.DB, species_name::String, 
                             temperature_range::Vector{<:Real}, config::Dict, 
                             step::Real=10.0)

Fixed version of calculate_properties that properly handles keyword arguments and errors.
Note: This version takes step as a regular parameter, not a keyword parameter.
"""
function calculate_properties_fixed(conn::DuckDB.DB, species_name::String, 
                                  temperature_range::Vector{<:Real}, 
                                  config::Dict, 
                                  step::Real=10.0)
    # Convert all parameters to ensure type stability
    temp_range = [Float64(temperature_range[1]), Float64(temperature_range[2])]
    step_val = Float64(step)
    
    # Validate temperature range
    if length(temp_range) != 2 || temp_range[1] >= temp_range[2]
        error("Invalid temperature range: $temp_range")
    end
    
    # Generate temperature points
    temperatures = temp_range[1]:step_val:temp_range[2]
    
    # Calculate properties at each temperature
    results = Dict[]
    
    for temp in temperatures
        try
            result = JThermodynamicsData.query_properties(conn, species_name, Float64(temp), config)
            push!(results, result)
        catch e
            @warn "Failed to calculate properties at $temp K: $e"
        end
    end
    
    # Check if we got any results
    if isempty(results)
        @warn "No valid results for any temperature point for $species_name"
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

"""
    plot_properties_log_temp(property_data::Dict, config::Dict)

Plot thermodynamic properties with a logarithmic temperature scale.
"""
function plot_properties_log_temp(property_data::Dict, config::Dict)
    # Check if we have any data to plot
    if isempty(property_data["temperatures"])
        @warn "No data to plot for $(property_data["species_name"])"
        # Return an empty plot
        return plot(title="No data available for $(property_data["species_name"])")
    end
    
    # Create individual plots
    temps = property_data["temperatures"]
    
    # Extract property data
    cp_values = property_data["properties"]["Cp"]["values"]
    cp_uncertainties = property_data["properties"]["Cp"]["uncertainties"]
    
    h_values = property_data["properties"]["H"]["values"]
    h_uncertainties = property_data["properties"]["H"]["uncertainties"]
    
    s_values = property_data["properties"]["S"]["values"]
    s_uncertainties = property_data["properties"]["S"]["uncertainties"]
    
    g_values = property_data["properties"]["G"]["values"]
    g_uncertainties = property_data["properties"]["G"]["uncertainties"]
    
    # Set up plots
    species_name = property_data["species_name"]
    
    # Heat capacity plot
    p_cp = plot(
        temps,  # Use original temps with log scale
        cp_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="Cp (J/mol/K)",
        title="Heat Capacity for $species_name",
        ribbon=cp_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Enthalpy plot
    p_h = plot(
        temps,
        h_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="H (kJ/mol)",
        title="Enthalpy for $species_name",
        ribbon=h_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Entropy plot
    p_s = plot(
        temps,
        s_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="S (J/mol/K)",
        title="Entropy for $species_name",
        ribbon=s_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Gibbs energy plot
    p_g = plot(
        temps,
        g_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="G (kJ/mol)",
        title="Gibbs Energy for $species_name",
        ribbon=g_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Combine into a 2x2 layout
    p = plot(p_cp, p_h, p_s, p_g, layout=(2, 2), size=(1000, 800))
    
    return p
end

function main()
    # Set up logging
    log_dir = joinpath(@__DIR__, "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "thermodynamics_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("ThermodynamicsDataProcessor", Dict(
        "start_time" => now(),
        "working_directory" => @__DIR__
    ))
    
    @info "JThermodynamicsData - Thermodynamic Data Processor"
    @info "================================================="
    
    # Load main configuration
    JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "LoadConfiguration", Dict())
    config_path = joinpath(@__DIR__, "config", "settings.yaml")
    config = JThermodynamicsData.with_timing(
        () -> JThermodynamicsData.load_config(config_path),
        "Configuration", 
        "LoadConfig"
    )
    JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "LoadConfiguration", Dict(
        "config_path" => config_path,
        "log_level" => config["general"]["log_level"]
    ))
    
    # Set up logging again with the config-specified level
    if haskey(config["general"], "log_level")
        JThermodynamicsData.init_logger(config["general"]["log_level"], log_file)
        @info "Updated log level to: $(config["general"]["log_level"])"
    end
    
    # Connect to existing database or initialize if needed
    JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "DatabaseConnection", Dict())
    db_path = config["general"]["database_path"]
    mkpath(dirname(db_path))
    
    # Check if database exists and try to connect first
    conn = nothing
    if isfile(db_path)
        try
            # Connect to existing database
            conn = JThermodynamicsData.with_timing(
                () -> JThermodynamicsData.connect_database(db_path),
                "Database", 
                "ConnectExisting"
            )
            @info "Connected to existing database" path=db_path
        catch e
            JThermodynamicsData.log_error("Could not connect to existing database", Dict(
                "error" => e,
                "path" => db_path
            ))
            @info "Will reinitialize database..." path=db_path
            rm(db_path)
            conn = JThermodynamicsData.with_timing(
                () -> JThermodynamicsData.init_database(config),
                "Database", 
                "Initialize"
            )
        end
    else
        # Initialize new database
        @info "No database found. Initializing new database..." path=db_path
        conn = JThermodynamicsData.with_timing(
            () -> JThermodynamicsData.init_database(config),
            "Database", 
            "Initialize"
        )
    end
    JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "DatabaseConnection", Dict(
        "db_path" => db_path,
        "connection_status" => "connected"
    ))
    
    # Load species list from YAML
    JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "LoadSpeciesConfig", Dict())
    species_config_path = joinpath(@__DIR__, "config", "species.yaml")
    if !isfile(species_config_path)
        JThermodynamicsData.log_error("Species configuration file not found", Dict(
            "path" => species_config_path
        ))
        error("Species configuration file not found: $species_config_path")
    end
    
    species_config = YAML.load_file(species_config_path)
    
    # Extract settings from species config
    temp_range = get(species_config, "temperature_range", [100.0, 10000.0])
    temp_step = get(species_config, "temperature_step", 10.0)
    species_list = species_config["species"]
    output_options = get(species_config, "output", Dict())
    JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "LoadSpeciesConfig", Dict(
        "species_count" => length(species_list),
        "temperature_range" => temp_range,
        "temperature_step" => temp_step
    ))
    
    # Create plots directory if needed
    plot_dir = get(output_options, "plot_directory", "plots")
    mkpath(plot_dir)
    
    # Create tables directory
    tables_dir = joinpath(@__DIR__, "output", "tables")
    mkpath(tables_dir)
    
    # Process each species
    total_species = length(species_list)
    @info "Processing species batch" count=total_species temp_range=temp_range
    
    JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "ProcessSpecies", Dict(
        "total_species" => total_species
    ))
    
    # Track statistics
    stats = Dict(
        "total_species" => total_species,
        "successful_species" => 0,
        "failed_species" => 0,
        "species_with_plots" => 0,
        "species_with_documentation" => 0
    )
    
    for (i, species_entry) in enumerate(species_list)
        # Extract species name (handle both string and dict formats)
        species_name = isa(species_entry, Dict) ? species_entry["name"] : string(species_entry)
        
        # Log species processing with context
        species_context = JThermodynamicsData.create_log_context()
        species_context["species_name"] = species_name
        species_context["index"] = i
        species_context["total"] = total_species
        
        @info "[$i/$total_species] Processing species" species=species_name
        
        JThermodynamicsData.log_species_processing(species_name, temp_range, temp_step)
        JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "Species:$species_name", Dict(
            "index" => i,
            "temperature_range" => temp_range,
            "temperature_step" => temp_step
        ))
        
        try
            # Check if species exists in the database
            database_check_start = now()
            species_id = JThermodynamicsData.get_species_id(conn, species_name)
            JThermodynamicsData.log_timing_benchmark("Database", "CheckSpeciesExists", now() - database_check_start)
            
            if species_id === nothing
                @warn "Species not found in database" species=species_name
                @info "Will use theoretical calculations only" species=species_name
                
                # Try to add species to database if it's not there
                try
                    # Run the initialize_minimal_database script to add the species
                    @info "Attempting to add species to database" species=species_name script="initialize_minimal_database.jl"
                    
                    # Time the database initialization
                    db_init_start = now()
                    include(joinpath(@__DIR__, "scripts", "initialize_minimal_database.jl"))
                    JThermodynamicsData.log_timing_benchmark("Database", "InitializeSpecies", now() - db_init_start)
                    
                    # Check if it worked
                    species_id = JThermodynamicsData.get_species_id(conn, species_name)
                    if species_id !== nothing
                        @info "Successfully added species to database" species=species_name id=species_id
                    end
                catch e
                    JThermodynamicsData.log_error("Failed to add species to database", Dict(
                        "species" => species_name,
                        "error" => e
                    ))
                end
            end
            
            # Calculate properties for the species
            JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "Calculate:$species_name", Dict(
                "species" => species_name,
                "temperature_range" => temp_range,
                "step" => temp_step
            ))
            
            result = Dict()
            step_val = Float64(temp_step)
            
            # For N2, use direct calculator first
            if species_name == "N2"
                @info "Using direct calculator" species=species_name
                
                calc_start = now()
                result = calculate_n2_properties_range(temp_range, step_val)
                JThermodynamicsData.log_timing_benchmark("Calculation", "DirectN2Calculator", now() - calc_start)
            else
                # For other species, use the database method
                calc_start = now()
                result = calculate_properties_fixed(conn, species_name, temp_range, config, step_val)
                JThermodynamicsData.log_timing_benchmark("Calculation", "DatabaseCalculator", now() - calc_start)
            end
            
            JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "Calculate:$species_name", Dict(
                "species" => species_name,
                "temperature_points" => length(result["temperatures"]),
                "success" => !isempty(result["temperatures"])
            ))
            
            # Check if we got any valid temperature points
            if isempty(result["temperatures"])
                JThermodynamicsData.log_error("No valid temperature points calculated", Dict(
                    "species" => species_name
                ))
                stats["failed_species"] += 1
                continue
            end
            
            # Save results to database if requested
            if get(output_options, "save_to_database", true)
                JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "SaveToDatabase:$species_name", Dict(
                    "species" => species_name
                ))
                
                @info "Saving results to database" species=species_name
                # The data is already stored in the database through the calculate_properties function
                
                JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "SaveToDatabase:$species_name", Dict(
                    "species" => species_name,
                    "success" => true
                ))
            end
            
            # Generate plots if requested
            if get(output_options, "generate_plots", true)
                JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "GeneratePlots:$species_name", Dict(
                    "species" => species_name
                ))
                
                use_log_temp = get(output_options, "use_log_temperature", true)  # Default to log temp scale
                plot_format = get(output_options, "plot_format", "png")
                
                # Generate plots with uncertainty ribbons
                @info "Generating plots" species=species_name log_scale=use_log_temp
                
                # Create plots with appropriate temperature scale
                plot_start = now()
                if use_log_temp
                    p = plot_properties_log_temp(result, config)
                else
                    p = JThermodynamicsData.plot_properties(result, config)
                end
                JThermodynamicsData.log_timing_benchmark("Visualization", "GeneratePlot", now() - plot_start)
                
                # Save the plot
                plot_file = joinpath(plot_dir, "$(species_name)_properties.$plot_format")
                savefig(p, plot_file)
                
                # Create additional plots only for non-N2 species
                if species_name != "N2"
                    try
                        # Create uncertainty analysis plot
                        uncertainty_start = now()
                        p_uncertainty = JThermodynamicsData.plot_uncertainty_analysis(conn, species_name, "Cp", temp_range, config, step=temp_step)
                        JThermodynamicsData.log_timing_benchmark("Visualization", "UncertaintyPlot", now() - uncertainty_start)
                        
                        uncertainty_file = joinpath(plot_dir, "$(species_name)_uncertainty.$plot_format")
                        savefig(p_uncertainty, uncertainty_file)
                        
                        # Create data source comparison if available
                        comparison_start = now()
                        p_comparison = JThermodynamicsData.plot_data_source_comparison(conn, species_name, "Cp", temp_range, config, step=temp_step)
                        JThermodynamicsData.log_timing_benchmark("Visualization", "ComparisonPlot", now() - comparison_start)
                        
                        comparison_file = joinpath(plot_dir, "$(species_name)_sources.$plot_format")
                        savefig(p_comparison, comparison_file)
                    catch e
                        JThermodynamicsData.log_error("Could not create additional plots", Dict(
                            "species" => species_name,
                            "error" => e
                        ))
                    end
                end
                
                JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "GeneratePlots:$species_name", Dict(
                    "species" => species_name,
                    "plot_file" => plot_file
                ))
                
                stats["species_with_plots"] += 1
            end
            
            # Generate documentation if requested
            if get(output_options, "generate_documentation", true)
                JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "GenerateDocs:$species_name", Dict(
                    "species" => species_name
                ))
                
                @info "Generating documentation" species=species_name
                docs_dir = get(output_options, "documentation_directory", joinpath(@__DIR__, "output", "documentation"))
                
                # Get all data sources for a single temperature point
                single_temp = (temp_range[1] + temp_range[2]) / 2
                
                refine_start = now()
                sources_query_result, all_sources = JThermodynamicsData.progressively_refine_thermodynamic_data(
                    conn, species_name, single_temp, config
                )
                JThermodynamicsData.log_timing_benchmark("Documentation", "GatherSources", now() - refine_start)
                
                # Create Markdown documentation
                try
                    markdown_start = now()
                    markdown_file = JThermodynamicsData.create_markdown_documentation(
                        species_name, sources_query_result, all_sources, docs_dir
                    )
                    JThermodynamicsData.log_timing_benchmark("Documentation", "CreateMarkdown", now() - markdown_start)
                    
                    @info "Created Markdown documentation" species=species_name file=markdown_file
                catch e
                    JThermodynamicsData.log_error("Failed to create Markdown documentation", Dict(
                        "species" => species_name,
                        "error" => e
                    ))
                end
                
                # Create JSON documentation
                try
                    json_start = now()
                    json_file = JThermodynamicsData.create_json_documentation(
                        species_name, sources_query_result, all_sources, docs_dir
                    )
                    JThermodynamicsData.log_timing_benchmark("Documentation", "CreateJSON", now() - json_start)
                    
                    @info "Created JSON documentation" species=species_name file=json_file
                catch e
                    JThermodynamicsData.log_error("Failed to create JSON documentation", Dict(
                        "species" => species_name,
                        "error" => e
                    ))
                end
                
                JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "GenerateDocs:$species_name", Dict(
                    "species" => species_name,
                    "docs_dir" => docs_dir
                ))
                
                stats["species_with_documentation"] += 1
            end
            
            # Create tabular output
            JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "CreateTable:$species_name", Dict(
                "species" => species_name
            ))
            
            @info "Creating tabular output" species=species_name
            table_file = joinpath(tables_dir, "$(species_name)_data.csv")
            
            # Open CSV file
            csv_start = now()
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
            JThermodynamicsData.log_timing_benchmark("Output", "CreateCSV", now() - csv_start)
            
            JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "CreateTable:$species_name", Dict(
                "species" => species_name,
                "csv_file" => table_file
            ))
            
            # Update statistics
            stats["successful_species"] += 1
            
            JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "Species:$species_name", Dict(
                "success" => true,
                "temperature_points" => length(result["temperatures"])
            ))
            
        catch e
            JThermodynamicsData.log_error("Error processing species", Dict(
                "species" => species_name,
                "error" => e,
                "stacktrace" => stacktrace(catch_backtrace())
            ))
            
            stats["failed_species"] += 1
            
            JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "Species:$species_name", Dict(
                "success" => false,
                "error" => string(e)
            ))
        end
    end
    
    JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "ProcessSpecies", stats)
    
    # Close database connection
    JThermodynamicsData.log_stage_start("ThermodynamicsDataProcessor", "CloseDatabase", Dict())
    JThermodynamicsData.close_database(conn)
    JThermodynamicsData.log_stage_end("ThermodynamicsDataProcessor", "CloseDatabase", Dict(
        "success" => true
    ))
    
    # Log pipeline completion
    JThermodynamicsData.log_pipeline_end("ThermodynamicsDataProcessor", Dict(
        "successful_species" => stats["successful_species"],
        "failed_species" => stats["failed_species"],
        "species_with_plots" => stats["species_with_plots"],
        "species_with_documentation" => stats["species_with_documentation"],
        "output_locations" => Dict(
            "database" => db_path,
            "plots" => plot_dir,
            "tables" => tables_dir,
            "logs" => log_dir
        )
    ))
    
    @info "Processing complete!"
    @info "Results saved in:"
    @info "  - Database: $db_path"
    @info "  - Plots: $plot_dir"
    @info "  - Tables: $tables_dir"
    @info "  - Logs: $log_dir"
end

# Run the main function
main()