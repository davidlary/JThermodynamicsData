#!/usr/bin/env julia

"""
Database Update Script

This script automatically updates the thermodynamic database from online sources.
It downloads the latest data from each configured source and updates the database.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using DuckDB
using DataFrames
using YAML
using Dates
using Printf

function main()
    println("JThermodynamicsData Database Update")
    println("==================================")
    
    # Load configuration
    config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
    config = load_config(config_path)
    
    # Initialize database
    db_path = config["general"]["database_path"]
    conn = init_database(config)
    
    # Create cache directory
    cache_dir = config["general"]["cache_directory"]
    mkpath(cache_dir)
    
    # Download and update data from each source
    update_all_sources(conn, config)
    
    # Perform database optimization
    optimize_database(conn)
    
    # Close database connection
    close_database(conn)
    
    println("\nDatabase update completed!")
    println("Last update: $(now())")
end

function update_all_sources(conn::DuckDB.DB, config::Dict)
    # Get data sources from configuration
    data_sources = config["data_sources"]
    
    # Sort by priority (lowest to highest)
    sort!(data_sources, by = source -> source["priority"])
    
    # Fetch and update each source
    for source in data_sources
        if source["enabled"]
            update_source(conn, source, config)
        else
            println("Skipping disabled source: $(source["name"])")
        end
    end
end

function update_source(conn::DuckDB.DB, source::Dict, config::Dict)
    name = source["name"]
    format = source["format"]
    reliability_score = get(source, "reliability_score", 0.0)
    
    println("\nUpdating source: $(name) (format: $(format))")
    
    # Create source-specific cache directory
    cache_dir = joinpath(config["general"]["cache_directory"], lowercase(name))
    mkpath(cache_dir)
    
    # Check if we need to update this source
    if should_update_source(conn, source)
        try
            # Fetch data for this source
            data_path = fetch_data_source(source, config["general"]["cache_directory"])
            
            if isempty(data_path) || !isfile(data_path)
                println("  No data file found for source: $(name)")
                return
            end
            
            # Process data and update database based on format
            species_count = update_database_from_source(conn, data_path, name, format, reliability_score)
            
            # Update last updated timestamp
            update_source_timestamp(conn, name)
            
            println("  Updated $(species_count) species for source: $(name)")
            
        catch e
            println("  Error updating source $(name): $(e)")
        end
    else
        println("  Source is up to date: $(name)")
    end
end

function should_update_source(conn::DuckDB.DB, source::Dict)
    name = source["name"]
    refresh_interval = get(source, "refresh_interval_days", 90) * 24 * 60 * 60  # Convert days to seconds
    
    # Query the last update time for this source
    query = """
    SELECT last_updated 
    FROM data_sources 
    WHERE name = ?
    """
    
    result = DuckDB.execute(conn, query, [name])
    df = DataFrame(result)
    
    if size(df, 1) > 0 && !ismissing(df[1, :last_updated])
        last_updated = df[1, :last_updated]
        
        # Calculate age in seconds
        age = time() - Dates.datetime2unix(last_updated)
        
        # Check if it's time to update
        if age < refresh_interval
            days_ago = round(age / (24 * 60 * 60), digits=1)
            println("  Last updated $(days_ago) days ago (interval: $(refresh_interval / (24 * 60 * 60)) days)")
            return false
        else
            days_ago = round(age / (24 * 60 * 60), digits=1)
            println("  Last updated $(days_ago) days ago - needs update")
            return true
        end
    end
    
    # If no record exists or last_updated is NULL, update is needed
    println("  No previous update record found - update needed")
    return true
end

function update_source_timestamp(conn::DuckDB.DB, source_name::String)
    # Update the last_updated timestamp for the source
    query = """
    UPDATE data_sources 
    SET last_updated = CURRENT_TIMESTAMP 
    WHERE name = ?
    """
    
    DuckDB.execute(conn, query, [source_name])
end

function update_database_from_source(conn::DuckDB.DB, data_path::String, name::String, format::String, reliability_score::Float64)
    # Process data based on format and update database
    species_count = 0
    
    if format == "nasa7" || format == "nasa7-9"
        species_count = nasa_to_database(conn, data_path, name, reliability_score)
    elseif format == "burcat"
        species_count = burcat_to_database(conn, data_path, name, reliability_score)
    elseif format == "chemkin"
        species_count = chemkin_to_database(conn, data_path, name, reliability_score)
    elseif format == "tabular" || format == "janaf"
        # For JANAF, we typically process one file per species
        # Here we assume data_path points to a directory with JANAF files
        if isdir(data_path)
            # Process each JANAF file in the directory
            for file in readdir(data_path, join=true)
                if endswith(lowercase(file), ".txt") || endswith(lowercase(file), ".dat")
                    try
                        species_name = janaf_to_database(conn, file, name, reliability_score)
                        if !isempty(species_name)
                            species_count += 1
                        end
                    catch e
                        println("  Error processing JANAF file $(file): $(e)")
                    end
                end
            end
        else
            # Single JANAF file
            species_name = janaf_to_database(conn, data_path, name, reliability_score)
            if !isempty(species_name)
                species_count += 1
            end
        end
    elseif format == "thermoml"
        species_count = thermoml_to_database(conn, data_path, name, reliability_score)
    elseif format == "tde"
        species_count = tde_to_database(conn, data_path, name, reliability_score)
    elseif format == "atct"
        species_count = atct_to_database(conn, data_path, name, reliability_score)
    else
        println("  Unsupported format: $(format)")
    end
    
    return species_count
end

function optimize_database(conn::DuckDB.DB)
    println("\nOptimizing database...")
    
    # Vacuum the database to reclaim space
    DuckDB.execute(conn, "VACUUM")
    
    # Analyze tables for query optimization
    DuckDB.execute(conn, "ANALYZE species")
    DuckDB.execute(conn, "ANALYZE thermodynamic_data")
    DuckDB.execute(conn, "ANALYZE data_sources")
    
    println("  Database optimization complete")
    
    # Print database statistics
    println("\nDatabase statistics:")
    
    # Count species
    result = DuckDB.execute(conn, "SELECT COUNT(*) as count FROM species")
    df = DataFrame(result)
    species_count = df[1, :count]
    println("  Total species: $(species_count)")
    
    # Count thermodynamic data entries
    result = DuckDB.execute(conn, "SELECT COUNT(*) as count FROM thermodynamic_data")
    df = DataFrame(result)
    data_count = df[1, :count]
    println("  Total data entries: $(data_count)")
    
    # Count data entries per source
    result = DuckDB.execute(conn, """
    SELECT td.data_source, COUNT(*) as count 
    FROM thermodynamic_data td 
    GROUP BY td.data_source
    ORDER BY count DESC
    """)
    df = DataFrame(result)
    
    println("  Data entries by source:")
    for row in eachrow(df)
        println("    $(row.data_source): $(row.count)")
    end
end

# Run the main function
main()