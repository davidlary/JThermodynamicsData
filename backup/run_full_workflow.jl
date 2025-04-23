#!/usr/bin/env julia

"""
Run Complete Thermodynamic Data Workflow

This script runs the complete thermodynamic data workflow:
1. Download all data sources from original locations 
2. Process and fetch data for all species from all sources
3. Generate plots for all species with hierarchical selection
4. Sync all data to DuckDB database

This gives you a complete setup with all data sources properly processed.
"""

using Pkg
Pkg.activate(@__DIR__)
import Dates

# Record the start time
start_time = Dates.now()

println("\n")
println("=" ^ 80)
println("JThermodynamicsData - Complete Workflow")
println("=" ^ 80)
println("Starting at: $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
println("\n")

# Step 1: Download ALL data sources
println("\n" * "=" ^ 60)
println("STEP 1: DOWNLOADING ALL DATA SOURCES")
println("=" ^ 60)

download_script = joinpath(@__DIR__, "scripts", "download_all_sources.jl")
if !isfile(download_script)
    error("Download script not found: $download_script")
end

println("Running: $download_script")
include(download_script)

# Step 2: Process and fetch data for all species from all sources
println("\n" * "=" ^ 60)
println("STEP 2: PROCESSING ALL DATA SOURCES")
println("=" ^ 60)

fetch_script = joinpath(@__DIR__, "scripts", "fetch_all_sources.jl")
if !isfile(fetch_script)
    error("Fetch script not found: $fetch_script")
end

println("Running: $fetch_script")
include(fetch_script)

# Step 3: Generate plots for all species with hierarchical selection
println("\n" * "=" ^ 60)
println("STEP 3: GENERATING PLOTS FOR ALL SPECIES")
println("=" ^ 60)

plot_script = joinpath(@__DIR__, "run_all_species_plots.jl")
if !isfile(plot_script)
    error("Plot script not found: $plot_script")
end

println("Running: $plot_script")
include(plot_script)

# Step 4: Sync all data to DuckDB database
println("\n" * "=" ^ 60)
println("STEP 4: SYNCING TO DATABASE")
println("=" ^ 60)

sync_script = joinpath(@__DIR__, "scripts", "sync_json_to_database.jl")
if !isfile(sync_script)
    error("Sync script not found: $sync_script")
end

println("Running: $sync_script")
include(sync_script)

# Record the end time and calculate duration
end_time = Dates.now()
duration = Dates.canonicalize(Dates.CompoundPeriod(end_time - start_time))

println("\n")
println("=" ^ 80)
println("WORKFLOW COMPLETE")
println("=" ^ 80)
println("Started at: $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
println("Finished at: $(Dates.format(end_time, "yyyy-mm-dd HH:MM:SS"))")
println("Total duration: $duration")
println("\n")
println("ALL thermodynamic data has been:")
println("1. Downloaded from original sources")
println("2. Processed for ALL species")
println("3. Visualized with hierarchical source selection")
println("4. Saved to DuckDB database")
println("\n")
println("Outputs can be found in:")
println("- Plots: $(joinpath(@__DIR__, "plots"))")
println("- Data: $(joinpath(@__DIR__, "data", "species"))")
println("- Database: $(joinpath(@__DIR__, "data", "thermodynamics.duckdb"))")
println("- Documentation: $(joinpath(@__DIR__, "output", "docs"))")
println("=" ^ 80)