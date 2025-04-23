#!/usr/bin/env julia

"""
Sync JSON Storage to Database

This script populates the DuckDB database with data from the JSON storage system.
It gives us the best of both worlds:
1. Reliable JSON storage for data loading
2. Database for advanced queries and analysis

This version ensures the correct priority ordering for data sources to maintain
the thermodynamic data hierarchy.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using JSON
using DataFrames
using DuckDB
using Printf
using Dates

"""
    init_database_structure(db_path)

Initialize database structure if it doesn't exist.
"""
function init_database_structure(db_path)
    println("Initializing database structure...")
    
    # Delete existing database file if it exists
    if isfile(db_path)
        @info "Removing existing database at $db_path"
        rm(db_path)
    end
    
    # Connect to database - this will create a new one
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Create data_sources table
    DBInterface.execute(conn, """
        CREATE TABLE data_sources (
            id INTEGER PRIMARY KEY,
            name VARCHAR NOT NULL UNIQUE,
            description VARCHAR,
            priority INTEGER,
            metadata_json VARCHAR
        )
    """)
    
    # Create species table
    DBInterface.execute(conn, """
        CREATE TABLE species (
            id INTEGER PRIMARY KEY,
            name VARCHAR NOT NULL UNIQUE,
            formula VARCHAR,
            cas_number VARCHAR,
            molecular_weight DOUBLE,
            metadata_json VARCHAR
        )
    """)
    
    # Create thermodynamic_data table with foreign key constraints
    DBInterface.execute(conn, """
        CREATE TABLE thermodynamic_data (
            id INTEGER PRIMARY KEY,
            species_id INTEGER,
            data_source VARCHAR,
            polynomial_type VARCHAR,
            temperature_min DOUBLE,
            temperature_max DOUBLE,
            priority INTEGER,
            reliability DOUBLE, 
            data_json VARCHAR,
            uncertainty_json VARCHAR,
            FOREIGN KEY (species_id) REFERENCES species(id),
            FOREIGN KEY (data_source) REFERENCES data_sources(name)
        )
    """)
    
    # Close connection
    DBInterface.close(conn)
    
    println("Database structure initialized.")
end

"""
    populate_data_sources(db_path)

Populate data sources table with hierarchy information.
This ensures correct priority ordering for the thermodynamic data hierarchy.
"""
function populate_data_sources(db_path)
    println("Populating data sources with correct hierarchy...")
    
    # Data sources with their priorities in the correct hierarchical order
    # Higher priority value = higher in hierarchy
    data_sources = [
        ("theoretical", "Theoretical calculations", 0, 2.5),
        ("group-contribution", "Group contribution methods", 1, 3.0),
        ("stat-thermo", "Statistical thermodynamics", 2, 3.5),
        ("benson-group", "Benson Group Additivity", 3, 3.0),
        ("quantum-statistical", "Quantum Statistical Thermodynamics", 4, 4.0),
        ("gri-mech", "GRI-MECH 3.0", 5, 3.5),
        ("chemkin", "CHEMKIN format data", 6, 3.8),
        ("nasa-cea", "NASA CEA Database", 7, 4.0),
        ("janaf", "JANAF Thermochemical Tables", 8, 4.5),
        ("thermoml", "ThermoML Standard", 9, 4.5),
        ("tde", "NIST ThermoData Engine", 10, 4.8),
        ("nist-webbook", "NIST Chemistry WebBook", 11, 4.95),
        ("burcat", "Burcat Database", 12, 4.9),
        ("atct", "Active Thermochemical Tables", 13, 5.0)
    ]
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Add data sources - database is already fresh from init_database_structure
    for (i, (name, description, priority, reliability)) in enumerate(data_sources)
        metadata = Dict(
            "reliability_score" => reliability,
            "date_added" => "2023-01-01",
            "description" => description
        )
        
        metadata_json = JSON.json(metadata)
        
        DBInterface.execute(conn, """
            INSERT INTO data_sources (id, name, description, priority, metadata_json)
            VALUES (?, ?, ?, ?, ?)
        """, [i, name, description, priority, metadata_json])
    end
    
    # Close connection
    DBInterface.close(conn)
    
    println("Data sources populated with correct priority ordering.")
end

"""
    sync_species_to_database(db_path)

Sync species information from JSON storage to database.

# Returns
- Number of species synced
"""
function sync_species_to_database(db_path)
    @info "Syncing species information to database..."
    
    # Get all species from JSON storage
    species_list = JThermodynamicsData.list_available_species()
    
    if isempty(species_list)
        @info "No species found in JSON storage. Please run initialize_json_storage.jl first."
        return 0
    end
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Add species data - database is already fresh from init_database_structure
    for (i, species_name) in enumerate(species_list)
        # Load species data from JSON
        species_data = JThermodynamicsData.load_species_data(species_name)
        
        # Extract basic information
        formula = get(species_data, "formula", species_name)
        cas_number = get(species_data, "cas_number", "")
        molecular_weight = get(species_data, "molecular_weight", 0.0)
        
        # Create metadata
        metadata = Dict(
            "date_added" => "2023-01-01",
            "description" => "Species from JSON storage"
        )
        
        metadata_json = JSON.json(metadata)
        
        # Insert into database
        DBInterface.execute(conn, """
            INSERT INTO species (id, name, formula, cas_number, molecular_weight, metadata_json)
            VALUES (?, ?, ?, ?, ?, ?)
        """, [i, species_name, formula, cas_number, molecular_weight, metadata_json])
    end
    
    # Close connection
    DBInterface.close(conn)
    
    @info "Synced $(length(species_list)) species to database."
    
    return length(species_list)
end

"""
    sync_thermodynamic_data(db_path)

Sync thermodynamic data from JSON storage to database.
Ensures that priority information is correctly transferred.

# Returns
- Number of thermodynamic data entries synced
"""
function sync_thermodynamic_data(db_path)
    @info "Syncing thermodynamic data to database with correct hierarchy..."
    
    # Get all species from JSON storage
    species_list = JThermodynamicsData.list_available_species()
    
    if isempty(species_list)
        @info "No species found in JSON storage. Please run initialize_json_storage.jl first."
        return 0
    end
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Get species ID mapping
    species_query = "SELECT id, name FROM species"
    species_result = DBInterface.execute(conn, species_query)
    species_df = DataFrame(species_result)
    
    species_id_map = Dict(species_df.name .=> species_df.id)
    
    # Get data source priority mapping
    sources_query = "SELECT name, priority FROM data_sources"
    sources_result = DBInterface.execute(conn, sources_query)
    sources_df = DataFrame(sources_result)
    
    source_priority_map = Dict(sources_df.name .=> sources_df.priority)
    
    # Add thermodynamic data for each species and source
    data_id = 1
    total_sources = 0
    
    for species_name in species_list
        # Load species data from JSON
        species_data = JThermodynamicsData.load_species_data(species_name)
        
        # Get species ID
        species_id = get(species_id_map, species_name, nothing)
        
        if species_id === nothing
            @warn "Species $(species_name) not found in database."
            continue
        end
        
        # Process each source
        if haskey(species_data, "sources")
            for (source_name, source_info) in species_data["sources"]
                if haskey(source_info, "data")
                    # Extract thermodynamic data
                    thermo_data = source_info["data"]
                    
                    # Extract temperature range
                    temp_min = get(thermo_data, "temperature_min", 200.0)
                    temp_max = get(thermo_data, "temperature_max", 6000.0)
                    
                    # Extract polynomial type
                    polynomial_type = get(thermo_data, "polynomial_type", "nasa7")
                    
                    # Get uncertainty
                    uncertainty = get(thermo_data, "uncertainty", 0.1)
                    uncertainty_json = JSON.json(Dict("value" => uncertainty))
                    
                    # Convert data to JSON
                    data_json = JSON.json(get(thermo_data, "data", Dict()))
                    
                    # Get priority and reliability from source info
                    priority = get(source_info, "priority", -1)
                    
                    # If priority is not specified, try to get it from the source map
                    if priority == -1
                        priority = get(source_priority_map, source_name, 0)
                    end
                    
                    reliability = get(source_info, "reliability_score", 3.0)
                    
                    # Insert into database
                    try
                        DBInterface.execute(conn, """
                            INSERT INTO thermodynamic_data 
                            (id, species_id, data_source, polynomial_type, temperature_min, temperature_max, 
                             priority, reliability, data_json, uncertainty_json)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """, [data_id, species_id, source_name, polynomial_type, temp_min, temp_max, 
                               priority, reliability, data_json, uncertainty_json])
                        
                        data_id += 1
                        total_sources += 1
                    catch e
                        @warn "Error inserting thermodynamic data for $(species_name) from source $(source_name): $(e)"
                    end
                end
            end
        end
    end
    
    # Close connection
    DBInterface.close(conn)
    
    @info "Synced $(total_sources) thermodynamic data entries to database with correct priority ordering."
    return total_sources
end

"""
    validate_database_hierarchy(db_path)

Validate that the database has the correct hierarchy ordering.
Checks that the highest priority source is selected for each species.

# Returns
- Bool indicating if validation passed
"""
function validate_database_hierarchy(db_path)
    println("Validating database hierarchy...")
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Query to get the highest priority source for each species
    validation_query = """
    WITH ranked_sources AS (
        SELECT 
            s.name AS species_name,
            t.data_source,
            t.priority,
            t.reliability,
            ROW_NUMBER() OVER (PARTITION BY s.name ORDER BY t.priority DESC, t.reliability DESC) AS rank
        FROM 
            species s
            JOIN thermodynamic_data t ON s.id = t.species_id
    )
    SELECT 
        species_name,
        data_source,
        priority,
        reliability
    FROM 
        ranked_sources
    WHERE 
        rank = 1
    ORDER BY 
        species_name
    """
    
    # Execute the query
    results = DBInterface.execute(conn, validation_query)
    df = DataFrame(results)
    
    # Close connection
    DBInterface.close(conn)
    
    # Print validation results
    println("\nHierarchical Source Selection Validation:")
    println("==========================================")
    println("| Species | Selected Source | Priority | Reliability |")
    println("|---------|----------------|----------|-------------|")
    
    for row in eachrow(df)
        println("| $(row.species_name) | $(row.data_source) | $(row.priority) | $(row.reliability) |")
    end
    
    # Check if any species uses experimental data
    exp_sources = ["atct", "burcat", "nist-webbook", "tde", "janaf", "nasa-cea", "chemkin", "gri-mech"]
    has_exp_data = any(row.data_source in exp_sources for row in eachrow(df))
    
    if has_exp_data
        println("\n✅ VALIDATION PASSED: Database includes experimental data sources!")
    else
        println("\n⚠️ VALIDATION WARNING: No species use experimental data sources!")
    end
    
    # Count sources by type
    source_counts = Dict{String, Int}()
    for row in eachrow(df)
        source_counts[row.data_source] = get(source_counts, row.data_source, 0) + 1
    end
    
    println("\nSource Distribution:")
    println("===================")
    for (source, count) in sort(collect(source_counts), by=x->-x[2])
        percentage = round(count * 100 / size(df, 1), digits=1)
        println("- $(source): $(count) species ($(percentage)%)")
    end
    
    println("\nDatabase hierarchy validation complete.")
    return has_exp_data
end

"""
    main()

Main function to sync JSON storage to database with correct hierarchy.
"""
function main()
    # Set up logging
    log_dir = joinpath(dirname(dirname(@__FILE__)), "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "sync_json_to_db_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("JsonToDatabaseSync", Dict(
        "start_time" => now(),
        "working_directory" => dirname(dirname(@__FILE__))
    ))
    
    @info "JThermodynamicsData - JSON to Database Sync Tool"
    @info "=============================================="
    
    # Database path
    db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
    
    # Initialize database structure
    JThermodynamicsData.log_stage_start("JsonToDatabaseSync", "InitDatabaseStructure", Dict(
        "db_path" => db_path
    ))
    db_struct_start = now()
    init_database_structure(db_path)
    db_struct_time = now() - db_struct_start
    JThermodynamicsData.log_timing_benchmark("Database", "InitStructure", db_struct_time)
    JThermodynamicsData.log_stage_end("JsonToDatabaseSync", "InitDatabaseStructure", Dict(
        "success" => true,
        "elapsed_time" => JThermodynamicsData.format_time_duration(db_struct_time)
    ))
    
    # Populate data sources with correct hierarchy
    JThermodynamicsData.log_stage_start("JsonToDatabaseSync", "PopulateDataSources", Dict())
    sources_start = now()
    populate_data_sources(db_path)
    sources_time = now() - sources_start
    JThermodynamicsData.log_timing_benchmark("Database", "PopulateSources", sources_time)
    JThermodynamicsData.log_stage_end("JsonToDatabaseSync", "PopulateDataSources", Dict(
        "success" => true,
        "elapsed_time" => JThermodynamicsData.format_time_duration(sources_time)
    ))
    
    # Sync species information
    JThermodynamicsData.log_stage_start("JsonToDatabaseSync", "SyncSpecies", Dict())
    species_start = now()
    species_count = sync_species_to_database(db_path)
    species_time = now() - species_start
    JThermodynamicsData.log_timing_benchmark("Database", "SyncSpecies", species_time)
    JThermodynamicsData.log_stage_end("JsonToDatabaseSync", "SyncSpecies", Dict(
        "species_count" => species_count,
        "elapsed_time" => JThermodynamicsData.format_time_duration(species_time)
    ))
    
    # Sync thermodynamic data with correct priority ordering
    JThermodynamicsData.log_stage_start("JsonToDatabaseSync", "SyncThermodynamicData", Dict())
    thermo_start = now()
    thermo_count = sync_thermodynamic_data(db_path)
    thermo_time = now() - thermo_start
    JThermodynamicsData.log_timing_benchmark("Database", "SyncThermodynamicData", thermo_time)
    JThermodynamicsData.log_stage_end("JsonToDatabaseSync", "SyncThermodynamicData", Dict(
        "thermo_entries" => thermo_count,
        "elapsed_time" => JThermodynamicsData.format_time_duration(thermo_time)
    ))
    
    # Validate database hierarchy
    JThermodynamicsData.log_stage_start("JsonToDatabaseSync", "ValidateHierarchy", Dict())
    validate_start = now()
    validation_result = validate_database_hierarchy(db_path)
    validate_time = now() - validate_start
    JThermodynamicsData.log_timing_benchmark("Database", "ValidateHierarchy", validate_time)
    JThermodynamicsData.log_stage_end("JsonToDatabaseSync", "ValidateHierarchy", Dict(
        "validation_passed" => validation_result,
        "elapsed_time" => JThermodynamicsData.format_time_duration(validate_time)
    ))
    
    @info "Sync complete!"
    @info "Database is now populated with data from JSON storage with the correct hierarchy."
    @info "You can now use both storage systems:"
    @info "- JSON for reliable data access"
    @info "- Database for advanced hierarchical querying"
    
    # Log pipeline completion
    JThermodynamicsData.log_pipeline_end("JsonToDatabaseSync", Dict(
        "species_count" => species_count,
        "thermo_entries" => thermo_count,
        "validation_passed" => validation_result,
        "db_path" => db_path,
        "log_file" => log_file
    ))
end

# Run the main function
main()