#!/usr/bin/env julia

"""
Sync JSON Storage to Database

This script populates the DuckDB database with data from the JSON storage system.
It gives us the best of both worlds:
1. Reliable JSON storage for data loading
2. Database for advanced queries and analysis
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using JSON
using DataFrames
using DuckDB
using Printf

"""
    init_database_structure(db_path)

Initialize database structure if it doesn't exist.
"""
function init_database_structure(db_path)
    println("Initializing database structure...")
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Create data_sources table
    DBInterface.execute(conn, """
        CREATE TABLE IF NOT EXISTS data_sources (
            id INTEGER PRIMARY KEY,
            name VARCHAR NOT NULL UNIQUE,
            description VARCHAR,
            priority INTEGER,
            metadata_json VARCHAR
        )
    """)
    
    # Create species table
    DBInterface.execute(conn, """
        CREATE TABLE IF NOT EXISTS species (
            id INTEGER PRIMARY KEY,
            name VARCHAR NOT NULL UNIQUE,
            formula VARCHAR,
            cas_number VARCHAR,
            molecular_weight DOUBLE,
            metadata_json VARCHAR
        )
    """)
    
    # Create thermodynamic_data table
    DBInterface.execute(conn, """
        CREATE TABLE IF NOT EXISTS thermodynamic_data (
            id INTEGER PRIMARY KEY,
            species_id INTEGER,
            data_source VARCHAR,
            polynomial_type VARCHAR,
            temperature_min DOUBLE,
            temperature_max DOUBLE,
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
"""
function populate_data_sources(db_path)
    println("Populating data sources...")
    
    # Data sources with their priorities
    data_sources = [
        ("theoretical", "Theoretical calculations", 0, 2.5),
        ("group-contribution", "Group contribution methods", 1, 3.0),
        ("stat-thermo", "Statistical thermodynamics", 2, 3.5),
        ("gri-mech", "GRI-MECH 3.0", 3, 3.5),
        ("chemkin", "CHEMKIN format data", 4, 3.8),
        ("nasa-cea", "NASA CEA Database", 5, 4.0),
        ("janaf", "JANAF Thermochemical Tables", 6, 4.5),
        ("thermoml", "ThermoML Standard", 7, 4.5),
        ("tde", "NIST ThermoData Engine", 8, 4.8),
        ("nist-webbook", "NIST Chemistry WebBook", 9, 4.95),
        ("burcat", "Burcat Database", 10, 4.9),
        ("atct", "Active Thermochemical Tables", 11, 5.0)
    ]
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Clear existing data sources
    DBInterface.execute(conn, "DELETE FROM data_sources")
    
    # Add data sources
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
    
    println("Data sources populated.")
end

"""
    sync_species_to_database(db_path)

Sync species information from JSON storage to database.
"""
function sync_species_to_database(db_path)
    println("Syncing species information to database...")
    
    # Get all species from JSON storage
    species_list = JThermodynamicsData.list_available_species()
    
    if isempty(species_list)
        println("No species found in JSON storage. Please run initialize_json_storage.jl first.")
        return
    end
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Clear existing species data
    DBInterface.execute(conn, "DELETE FROM species")
    
    # Add species data
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
    
    println("Synced $(length(species_list)) species to database.")
end

"""
    sync_thermodynamic_data(db_path)

Sync thermodynamic data from JSON storage to database.
"""
function sync_thermodynamic_data(db_path)
    println("Syncing thermodynamic data to database...")
    
    # Get all species from JSON storage
    species_list = JThermodynamicsData.list_available_species()
    
    if isempty(species_list)
        println("No species found in JSON storage. Please run initialize_json_storage.jl first.")
        return
    end
    
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Clear existing thermodynamic data
    DBInterface.execute(conn, "DELETE FROM thermodynamic_data")
    
    # Get species ID mapping
    species_query = "SELECT id, name FROM species"
    species_result = DBInterface.execute(conn, species_query)
    species_df = DataFrame(species_result)
    
    species_id_map = Dict(species_df.name .=> species_df.id)
    
    # Add thermodynamic data for each species and source
    data_id = 1
    total_sources = 0
    
    for species_name in species_list
        # Load species data from JSON
        species_data = JThermodynamicsData.load_species_data(species_name)
        
        # Get species ID
        species_id = get(species_id_map, species_name, nothing)
        
        if species_id === nothing
            println("  Warning: Species $(species_name) not found in database.")
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
                    
                    # Insert into database
                    DBInterface.execute(conn, """
                        INSERT INTO thermodynamic_data 
                        (id, species_id, data_source, polynomial_type, temperature_min, temperature_max, data_json, uncertainty_json)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, [data_id, species_id, source_name, polynomial_type, temp_min, temp_max, data_json, uncertainty_json])
                    
                    data_id += 1
                    total_sources += 1
                end
            end
        end
    end
    
    # Close connection
    DBInterface.close(conn)
    
    println("Synced $(total_sources) thermodynamic data entries to database.")
end

"""
    main()

Main function to sync JSON storage to database.
"""
function main()
    println("JThermodynamicsData - JSON to Database Sync Tool")
    println("===============================================")
    
    # Database path
    db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
    
    # Initialize database structure
    init_database_structure(db_path)
    
    # Populate data sources
    populate_data_sources(db_path)
    
    # Sync species information
    sync_species_to_database(db_path)
    
    # Sync thermodynamic data
    sync_thermodynamic_data(db_path)
    
    println("\nSync complete!")
    println("Database is now populated with data from JSON storage.")
    println("You can now use both storage systems:")
    println("- JSON for reliable data access")
    println("- Database for advanced querying")
end

# Run the main function
main()