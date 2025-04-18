#!/usr/bin/env julia

# Script to initialize the thermodynamic database
# This script downloads data from all configured sources and initializes the database

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using DataFrames

println("Initializing JThermodynamicsData Database")
println("==========================================")

# Load configuration
config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
config = load_config(config_path)

# Create cache directory
cache_dir = config["general"]["cache_directory"]
mkpath(cache_dir)

# Initialize database
println("Initializing database...")
db_path = config["general"]["database_path"]
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

# Train machine learning model if enabled
ml_enabled = false

for method in config["theoretical_calculation"]["methods"]
    if method["name"] == "machine_learning" && method["enabled"]
        ml_enabled = true
        break
    end
end

if ml_enabled
    println("Training machine learning model...")
    try
        train_ml_model(conn, config)
        println("  Machine learning model trained successfully")
    catch e
        @error "Failed to train machine learning model: $e"
    end
end

# Pre-cache common species if enabled
if config["performance"]["cache_frequent_species"]
    println("Pre-caching common species...")
    try
        cache_frequent_species(conn, config)
        println("  Common species pre-cached")
    catch e
        @error "Failed to pre-cache common species: $e"
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

println("Database initialization complete!")