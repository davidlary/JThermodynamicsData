#!/usr/bin/env julia

"""
Reinitialize JThermodynamicsData Database

This script removes the existing database and reinitializes it from scratch.
Use this when you need to rebuild the database completely.

Usage:
    julia scripts/reinitialize_database.jl
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using DataFrames

println("Reinitializing JThermodynamicsData Database")
println("===========================================")

# Load configuration
config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
config = load_config(config_path)

# Get database path
db_path = config["general"]["database_path"]
println("Database path: $db_path")

# Delete existing database if it exists
if isfile(db_path)
    println("Removing existing database...")
    rm(db_path)
    println("  Existing database removed.")
else
    println("No existing database found. Creating new database.")
end

# Create cache directory
cache_dir = config["general"]["cache_directory"]
mkpath(cache_dir)

# Initialize database
println("Initializing new database...")
mkpath(dirname(db_path))
conn = init_database(config)

# Fetch all data sources
println("Fetching data sources...")
source_paths = fetch_all_sources(config)

# Import data from each source
data_counts = Dict{String, Int}()

for (source_name, source_path) in source_paths
    println("Importing data from $source_name...")
    
    # Get source configuration
    source_config = nothing
    reliability_score = 0.0
    
    for source in config["data_sources"]
        if source["name"] == source_name
            source_config = source
            reliability_score = get(source, "reliability_score", 0.0)
            break
        end
    end
    
    if source_config === nothing
        @warn "Source configuration not found for $source_name"
        continue
    end
    
    # Import based on format
    format = source_config["format"]
    count = 0
    
    try
        if format == "chemkin"
            count = chemkin_to_database(conn, source_path, source_name, reliability_score)
        elseif format == "nasa7-9"
            # Implement import for NASA CEA format
            @warn "NASA CEA import not yet implemented"
        elseif format == "tabular" || format == "janaf"
            # Implement import for JANAF format
            @warn "JANAF import not yet implemented"
        elseif format == "thermoml"
            # Implement import for ThermoML format
            @warn "ThermoML import not yet implemented"
        elseif format == "tde"
            # Implement import for TDE format
            @warn "TDE import not yet implemented"
        elseif format == "burcat"
            count = burcat_to_database(conn, source_path, source_name, reliability_score)
        elseif format == "atct"
            # Implement import for ATcT format
            @warn "ATcT import not yet implemented"
        else
            @warn "Unknown format: $format"
        end
        
        data_counts[source_name] = count
        println("  Imported $count species from $source_name")
    catch e
        @error "Failed to import data from $source_name: $e"
    end
end

# Print summary
println("\nDatabase Initialization Summary")
println("------------------------------")
println("Database path: $db_path")
println("Data sources imported:")

for (source_name, count) in data_counts
    println("  $source_name: $count species")
end

# Count total species
total_species_query = "SELECT COUNT(DISTINCT id) AS count FROM species"
total_species_result = DuckDB.execute(conn, total_species_query)
total_species_df = DataFrame(total_species_result)
total_species = total_species_df[1, :count]

println("Total species in database: $total_species")

# Vacuum and close the database
println("Finalizing database...")
vacuum_database(conn)
close_database(conn)

println("Database reinitialization complete!")