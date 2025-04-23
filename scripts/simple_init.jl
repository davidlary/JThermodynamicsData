#!/usr/bin/env julia

# Simple initialization script for the thermodynamic database
# This script creates a minimal database for the hierarchical demo

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using DuckDB
using DataFrames
using JSON
using YAML

println("Simple Initialization of JThermodynamicsData Database")
println("====================================================")

# Initialize database
db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
mkpath(dirname(db_path))

# Make sure cache directories exist
data_dir = joinpath(dirname(dirname(@__FILE__)), "data")
mkpath(joinpath(data_dir, "cache", "chemkin"))
mkpath(joinpath(data_dir, "cache", "burcat"))
mkpath(joinpath(data_dir, "cache", "janaf"))
mkpath(joinpath(data_dir, "cache", "nasa-cea"))
mkpath(joinpath(data_dir, "cache", "thermoml"))
mkpath(joinpath(data_dir, "cache", "tde"))
mkpath(joinpath(data_dir, "cache", "atct"))
mkpath(joinpath(data_dir, "external", "demo"))

# Remove existing database if it exists
if isfile(db_path)
    println("Removing existing database at $db_path")
    rm(db_path)
end

println("Creating new database at $db_path")

# Initialize with default config
config, db = initialize()

# Create sequences for auto-increment if not already created
println("Creating sequences...")
conn = JThermodynamicsData.get_connection(db)

# Drop sequences if they exist and recreate them
println("  Resetting sequences...")
try
    DuckDB.execute(conn, "DROP SEQUENCE IF EXISTS species_id_seq")
    DuckDB.execute(conn, "DROP SEQUENCE IF EXISTS thermodynamic_data_id_seq")
    DuckDB.execute(conn, "DROP SEQUENCE IF EXISTS data_sources_id_seq")
    DuckDB.execute(conn, "DROP SEQUENCE IF EXISTS data_cache_id_seq")
    DuckDB.execute(conn, "DROP SEQUENCE IF EXISTS database_version_id_seq")
catch e
    println("  Warning: Failed to drop sequences: $e")
end

# Create fresh sequences
println("  Creating new sequences...")
DuckDB.execute(conn, "CREATE SEQUENCE species_id_seq")
DuckDB.execute(conn, "CREATE SEQUENCE thermodynamic_data_id_seq")
DuckDB.execute(conn, "CREATE SEQUENCE data_sources_id_seq")
DuckDB.execute(conn, "CREATE SEQUENCE data_cache_id_seq")
DuckDB.execute(conn, "CREATE SEQUENCE database_version_id_seq")

# Setup hierarchical sources
println("Setting up data sources for hierarchical selection...")
sources = [
    ("theoretical", "Theoretical calculations", 1, 3.0),
    ("gri-mech", "GRI-MECH 3.0", 3, 3.5),
    ("burcat", "Burcat database", 7, 4.2),
    ("janaf", "JANAF tables", 8, 4.5),
    ("atct", "Active Thermochemical Tables", 9, 5.0)
]

# Add each source
for (i, (name, description, priority, reliability)) in enumerate(sources)
    println("  Adding source: $name (priority: $priority)")
    
    # Create metadata with reliability score
    metadata = Dict(
        "format" => "nasa7",
        "enabled" => true,
        "reliability_score" => reliability,
        "refresh_interval_days" => 90
    )
    
    metadata_json = JSON.json(metadata)
    
    # Add source to database
    DuckDB.execute(conn, """
        INSERT INTO data_sources (id, name, description, priority, metadata_json)
        VALUES (nextval('data_sources_id_seq'), ?, ?, ?, ?)
    """, [name, description, priority, metadata_json])
end

# Import sample data - add some common species
println("Adding sample species data...")

# Add sample species
species_list = [
    ("N2", "N2", "7727-37-9", 28.0134),
    ("O2", "O2", "7782-44-7", 31.9988),
    ("H2O", "H2O", "7732-18-5", 18.0153),
    ("CO2", "CO2", "124-38-9", 44.01)
]

# Add each species directly through SQL
for (name, formula, cas, mw) in species_list
    println("  Adding species: $name")
    
    # Add species to the database
    species_query = """
    INSERT INTO species (id, name, formula, cas_number, molecular_weight)
    VALUES (nextval('species_id_seq'), ?, ?, ?, ?)
    RETURNING id
    """
    
    result = DuckDB.execute(conn, species_query, [name, formula, cas, mw])
    species_id = DataFrame(result)[1, :id]
    
    # For each source, add sample data at different reliability levels
    for (source_name, _, _, reliability) in sources
        # Create NASA-7 coefficients for this species with slightly different values per source
        # to demonstrate the hierarchical selection
        
        data_json = JSON.json(Dict(
            "low_temp" => Dict(
                "range" => [200.0, 1000.0],
                "coefficients" => [3.5, 0.001*reliability, 0.0, 0.0, 0.0, -1000.0, 10.0]
            ),
            "high_temp" => Dict(
                "range" => [1000.0, 6000.0],
                "coefficients" => [3.5, 0.001*reliability, 0.0, 0.0, 0.0, -1000.0, 10.0]
            )
        ))
        
        # Create uncertainty JSON
        uncertainty_json = JSON.json(Dict(
            "type" => "constant",
            "value" => 0.2 - reliability * 0.02  # Higher reliability = lower uncertainty
        ))
        
        # Add thermodynamic data for this source
        println("    Adding data from source: $source_name (reliability: $reliability)")
        
        data_query = """
        INSERT INTO thermodynamic_data 
        (id, species_id, data_source, polynomial_type, temperature_min, temperature_max, 
         reliability_score, data_json, uncertainty_json)
        VALUES (nextval('thermodynamic_data_id_seq'), ?, ?, ?, ?, ?, ?, ?, ?)
        """
        
        DuckDB.execute(conn, data_query, [
            species_id, 
            source_name, 
            "nasa7", 
            200.0, 
            6000.0, 
            reliability,
            data_json,
            uncertainty_json
        ])
    end
end

# Vacuum and close the database
println("Finalizing database...")
DuckDB.execute(conn, "VACUUM")
close_database(db)

println("Simple database initialization complete!")
println("Run the hierarchical demo with: julia run_hierarchical_demo.jl")