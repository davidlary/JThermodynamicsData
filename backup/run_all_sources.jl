#!/usr/bin/env julia

"""
Run All Data Sources - Complete Thermodynamic Pipeline

This script runs the complete thermodynamic data pipeline:
1. Download all data sources
2. Parse and process data
3. Generate plots for comparison
4. Sync to database for querying

This will ensure the most accurate source is used for each species.
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using Dates
using Printf

"""
    main()

Main function to run all steps of the thermodynamic data pipeline.
"""
function main()
    # Set up logging
    log_dir = joinpath(@__DIR__, "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "complete_pipeline_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("CompletePipeline", Dict(
        "start_time" => now(),
        "working_directory" => @__DIR__
    ))
    
    println("JThermodynamicsData - Complete Pipeline Runner")
    println("=============================================")
    
    # Step 1: Download all data sources
    println("\n== Step 1: Downloading All Data Sources ==")
    JThermodynamicsData.log_stage_start("CompletePipeline", "DownloadSources", Dict())
    
    download_start = now()
    include(joinpath(@__DIR__, "scripts", "download_all_sources.jl"))
    download_time = now() - download_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "DownloadSources", download_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "DownloadSources", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(download_time)
    ))
    
    # Step 2: Fetch and process data for all species
    println("\n== Step 2: Fetching Data from All Sources ==")
    JThermodynamicsData.log_stage_start("CompletePipeline", "FetchAllData", Dict())
    
    fetch_start = now()
    include(joinpath(@__DIR__, "scripts", "fetch_all_sources.jl"))
    fetch_time = now() - fetch_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "FetchAllData", fetch_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "FetchAllData", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(fetch_time)
    ))
    
    # Step 3: Generate plots for all species
    println("\n== Step 3: Generating Comparison Plots ==")
    JThermodynamicsData.log_stage_start("CompletePipeline", "GeneratePlots", Dict())
    
    plots_start = now()
    include(joinpath(@__DIR__, "run_all_species_plots.jl"))
    plots_time = now() - plots_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "GeneratePlots", plots_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "GeneratePlots", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(plots_time)
    ))
    
    # Step 4: Sync data to database
    println("\n== Step 4: Syncing to Database ==")
    JThermodynamicsData.log_stage_start("CompletePipeline", "SyncToDatabase", Dict())
    
    db_start = now()
    include(joinpath(@__DIR__, "scripts", "sync_json_to_database.jl"))
    db_time = now() - db_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "SyncToDatabase", db_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "SyncToDatabase", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(db_time)
    ))
    
    # Finish pipeline
    total_time = now() - download_start
    JThermodynamicsData.log_pipeline_end("CompletePipeline", Dict(
        "total_time" => JThermodynamicsData.format_time_duration(total_time),
        "data_dir" => joinpath(@__DIR__, "data"),
        "plots_dir" => joinpath(@__DIR__, "plots"),
        "db_path" => joinpath(@__DIR__, "data", "thermodynamics.duckdb"),
        "log_file" => log_file
    ))
    
    println("\n== Complete Pipeline Execution Summary ==")
    println("Total Execution Time: $(JThermodynamicsData.format_time_duration(total_time))")
    println("  - Download Sources: $(JThermodynamicsData.format_time_duration(download_time))")
    println("  - Fetch All Data: $(JThermodynamicsData.format_time_duration(fetch_time))")
    println("  - Generate Plots: $(JThermodynamicsData.format_time_duration(plots_time))")
    println("  - Sync to Database: $(JThermodynamicsData.format_time_duration(db_time))")
    println()
    println("All data has been processed, plotted, and synced to the database.")
    println("You can now use the run_thermodynamics.jl script to query properties.")
end

# Run the main function
main()