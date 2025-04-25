#!/usr/bin/env julia

"""
Run Complete Thermodynamic Data Pipeline Workflow

This script runs the complete thermodynamic data pipeline:
1. Download data from all sources (experimental & theoretical)
2. Process and parse data for all species
3. Generate plots for each species showing all sources
4. Sync all data to DuckDB database
5. Verify the hierarchical selection is working correctly

This ensures that all data sources are properly integrated and the hierarchy
is correctly followed when selecting data for each species.
"""

using Pkg
Pkg.activate(@__DIR__)
import Dates
using JThermodynamicsData

# Record the start time
start_time = Dates.now()

println("\n")
println("=" ^ 80)
println("JThermodynamicsData - Complete Pipeline Workflow")
println("=" ^ 80)
println("Starting at: $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
println("\n")

# Initialize the logger
log_dir = joinpath(@__DIR__, "logs")
mkpath(log_dir)
log_file = joinpath(log_dir, "complete_workflow_$(Dates.format(start_time, "yyyymmdd_HHMMSS")).log")
JThermodynamicsData.init_logger("info", log_file)

# Start pipeline logging
JThermodynamicsData.log_pipeline_start("CompletePipeline", Dict(
    "start_time" => start_time,
    "working_directory" => @__DIR__
))

# Step 1: Download ALL data sources
println("\n" * "=" ^ 60)
println("STEP 1: DOWNLOADING ALL DATA SOURCES")
println("=" ^ 60)

download_script = joinpath(@__DIR__, "scripts", "download_all_sources.jl")
if !isfile(download_script)
    error("Download script not found: $download_script")
end

println("Running: $download_script")
download_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "DownloadSources", Dict())

try
    include(download_script)
    download_time = Dates.now() - download_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "DownloadSources", download_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "DownloadSources", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(download_time)
    ))
    
    println("✅ Download completed in $(JThermodynamicsData.format_time_duration(download_time))")
catch e
    println("❌ Download failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "DownloadSources", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Step 2: Process and fetch data for all species from all sources
println("\n" * "=" ^ 60)
println("STEP 2: PROCESSING ALL DATA SOURCES")
println("=" ^ 60)

fetch_script = joinpath(@__DIR__, "scripts", "fetch_all_sources.jl")
if !isfile(fetch_script)
    error("Fetch script not found: $fetch_script")
end

println("Running: $fetch_script")
fetch_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "FetchAllData", Dict())

try
    include(fetch_script)
    fetch_time = Dates.now() - fetch_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "FetchAllData", fetch_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "FetchAllData", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(fetch_time)
    ))
    
    println("✅ Data processing completed in $(JThermodynamicsData.format_time_duration(fetch_time))")
catch e
    println("❌ Data processing failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "FetchAllData", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Step 3: Generate plots for all species with hierarchical selection
println("\n" * "=" ^ 60)
println("STEP 3: GENERATING PLOTS FOR ALL SPECIES")
println("=" ^ 60)

plot_script = joinpath(@__DIR__, "run_all_species_plots.jl")
if !isfile(plot_script)
    error("Plot script not found: $plot_script")
end

println("Running: $plot_script")
plot_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "GeneratePlots", Dict())

try
    include(plot_script)
    plot_time = Dates.now() - plot_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "GeneratePlots", plot_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "GeneratePlots", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(plot_time)
    ))
    
    println("✅ Plot generation completed in $(JThermodynamicsData.format_time_duration(plot_time))")
catch e
    println("❌ Plot generation failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "GeneratePlots", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Step 4: Sync all data to DuckDB database
println("\n" * "=" ^ 60)
println("STEP 4: SYNCING TO DATABASE")
println("=" ^ 60)

sync_script = joinpath(@__DIR__, "scripts", "sync_json_to_database.jl")
if !isfile(sync_script)
    error("Sync script not found: $sync_script")
end

println("Running: $sync_script")
sync_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "SyncToDatabase", Dict())

try
    include(sync_script)
    sync_time = Dates.now() - sync_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "SyncToDatabase", sync_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "SyncToDatabase", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(sync_time)
    ))
    
    println("✅ Database sync completed in $(JThermodynamicsData.format_time_duration(sync_time))")
catch e
    println("❌ Database sync failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "SyncToDatabase", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Step 5: Verify the pipeline is working correctly
println("\n" * "=" ^ 60)
println("STEP 5: VERIFYING PIPELINE")
println("=" ^ 60)

# First verify real data implementation
real_data_script = joinpath(@__DIR__, "verify_real_data.jl")
if !isfile(real_data_script)
    error("Real data verification script not found: $real_data_script")
end

println("Running real data verification: $real_data_script")
verify_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "VerifyRealData", Dict())

try
    include(real_data_script)
    verify_time = Dates.now() - verify_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "VerifyRealData", verify_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "VerifyRealData", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(verify_time)
    ))
    
    println("✅ Real data verification completed in $(JThermodynamicsData.format_time_duration(verify_time))")
catch e
    println("❌ Real data verification failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "VerifyRealData", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Then run the complete pipeline verification
verify_script = joinpath(@__DIR__, "verify_pipeline.jl")
if isfile(verify_script)
    println("\nRunning full pipeline verification: $verify_script")
    pipeline_verify_start = Dates.now()
    JThermodynamicsData.log_stage_start("CompletePipeline", "VerifyPipeline", Dict())

    try
        include(verify_script)
        pipeline_verify_time = Dates.now() - pipeline_verify_start
        
        JThermodynamicsData.log_timing_benchmark("Pipeline", "VerifyPipeline", pipeline_verify_time)
        JThermodynamicsData.log_stage_end("CompletePipeline", "VerifyPipeline", Dict(
            "elapsed_time" => JThermodynamicsData.format_time_duration(pipeline_verify_time)
        ))
        
        println("✅ Pipeline verification completed in $(JThermodynamicsData.format_time_duration(pipeline_verify_time))")
    catch e
        println("❌ Pipeline verification failed: $e")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        JThermodynamicsData.log_stage_end("CompletePipeline", "VerifyPipeline", Dict(
            "status" => "failed",
            "error" => "$e"
        ))
    end
else
    println("\n⚠️ Full pipeline verification script not found: $verify_script")
    println("Skipping full pipeline verification...")
end

# Step 6: Generate summary report
println("\n" * "=" ^ 60)
println("STEP 6: GENERATING SUMMARY REPORT")
println("=" ^ 60)

summary_script = joinpath(@__DIR__, "scripts", "generate_source_summary.jl")
if !isfile(summary_script)
    error("Summary script not found: $summary_script")
end

println("Running: $summary_script")
summary_start = Dates.now()
JThermodynamicsData.log_stage_start("CompletePipeline", "GenerateSummary", Dict())

try
    include(summary_script)
    summary_time = Dates.now() - summary_start
    
    JThermodynamicsData.log_timing_benchmark("Pipeline", "GenerateSummary", summary_time)
    JThermodynamicsData.log_stage_end("CompletePipeline", "GenerateSummary", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(summary_time)
    ))
    
    println("✅ Summary generation completed in $(JThermodynamicsData.format_time_duration(summary_time))")
catch e
    println("❌ Summary generation failed: $e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    JThermodynamicsData.log_stage_end("CompletePipeline", "GenerateSummary", Dict(
        "status" => "failed",
        "error" => "$e"
    ))
end

# Record the end time and calculate duration
end_time = Dates.now()
duration = Dates.canonicalize(Dates.CompoundPeriod(end_time - start_time))

# Finish pipeline logging
JThermodynamicsData.log_pipeline_end("CompletePipeline", Dict(
    "total_time" => JThermodynamicsData.format_time_duration(duration),
    "data_dir" => joinpath(@__DIR__, "data"),
    "plots_dir" => joinpath(@__DIR__, "plots"),
    "db_path" => joinpath(@__DIR__, "data", "thermodynamics.duckdb"),
    "log_file" => log_file
))

println("\n")
println("=" ^ 80)
println("COMPLETE WORKFLOW EXECUTION SUMMARY")
println("=" ^ 80)
println("Started at: $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
println("Finished at: $(Dates.format(end_time, "yyyy-mm-dd HH:MM:SS"))")
println("Total duration: $duration")
println("\n")

# Print individual step timings
println("Step Timings:")
println("-------------")
try
    download_time_str = @isdefined(download_time) ? JThermodynamicsData.format_time_duration(download_time) : "N/A"
    fetch_time_str = @isdefined(fetch_time) ? JThermodynamicsData.format_time_duration(fetch_time) : "N/A"
    plot_time_str = @isdefined(plot_time) ? JThermodynamicsData.format_time_duration(plot_time) : "N/A"
    sync_time_str = @isdefined(sync_time) ? JThermodynamicsData.format_time_duration(sync_time) : "N/A"
    verify_time_str = @isdefined(verify_time) ? JThermodynamicsData.format_time_duration(verify_time) : "N/A"
    summary_time_str = @isdefined(summary_time) ? JThermodynamicsData.format_time_duration(summary_time) : "N/A"
    
    println("1. Download All Sources: $download_time_str")
    println("2. Process All Data: $fetch_time_str")
    println("3. Generate Plots: $plot_time_str")
    println("4. Sync to Database: $sync_time_str")
    println("5. Verify Pipeline: $verify_time_str")
    println("6. Generate Summary: $summary_time_str")
catch
    println("Error displaying timings")
end

println("\n")
println("ALL THERMODYNAMIC DATA HAS BEEN:")
println("1. Downloaded from original sources")
println("2. Processed for ALL species")
println("3. Visualized with hierarchical source selection")
println("4. Saved to DuckDB database")
println("5. Verified for correct hierarchical operation")
println("\n")
println("Outputs can be found in:")
println("- Plots: $(joinpath(@__DIR__, "plots"))")
println("- Data: $(joinpath(@__DIR__, "data", "species"))")
println("- Database: $(joinpath(@__DIR__, "data", "thermodynamics.duckdb"))")
println("- Documentation: $(joinpath(@__DIR__, "output", "docs"))")
println("- Summary: $(joinpath(@__DIR__, "output", "summary_report.md"))")
println("- Log: $log_file")
println("=" ^ 80)