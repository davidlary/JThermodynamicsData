#!/usr/bin/env julia

# Initialize a minimal database with just N2 for debugging purposes

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using JSON
using DuckDB
using DataFrames

println("Initializing Minimal Database for Debugging")
println("===========================================")

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

# Add N2 to the database manually
println("Adding N2 to database...")

# Species table entry
species_id = 1
species_name = "N2"
species_formula = "N2"
species_cas = "7727-37-9"
molecular_weight = 28.0134

# Create molecular data for N2
metadata = Dict(
    "molecular_data" => Dict(
        "molecular_weight" => molecular_weight,
        "symmetry_number" => 2,
        "spin_multiplicity" => 1,
        "rotational_constants" => [1.998],  # cm^-1
        "vibrational_frequencies" => [2358.0]  # cm^-1
    )
)
metadata_json = JSON.json(metadata)

# First check if N2 already exists
check_query = "SELECT id FROM species WHERE name = ?"
check_result = DuckDB.execute(conn, check_query, [species_name])
check_df = DataFrame(check_result)

if size(check_df, 1) > 0
    println("N2 already exists in database with ID: $(check_df[1, :id])")
    species_id = check_df[1, :id]
    
    # Update the existing record with metadata
    update_query = """
    UPDATE species 
    SET metadata_json = ?
    WHERE id = ?
    """
    
    DuckDB.execute(conn, update_query, [metadata_json, species_id])
    println("Updated existing N2 record with metadata")
else
    # Add to species table - explicitly set the ID since it's not auto-incrementing
    species_query = """
    INSERT INTO species (id, name, formula, cas_number, molecular_weight, metadata_json)
    VALUES (?, ?, ?, ?, ?, ?)
    """
    
    DuckDB.execute(conn, species_query, [
        species_id, species_name, species_formula, species_cas, molecular_weight, metadata_json
    ])
    println("Inserted N2 with ID: $species_id")
end

# Add basic NASA-7 polynomial data for N2
# Using NASA-7 coefficients from GRI-Mech 3.0
nasa_coeff = [
    # Low temperature range (200-1000K)
    [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372],
    # High temperature range (1000-6000K)
    [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]
]

# JSON representation of coefficients
data_json = JSON.json(Dict(
    "coefficients" => nasa_coeff,
    "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]]
))

# Check if data already exists for N2
data_check_query = "SELECT id FROM thermodynamic_data WHERE species_id = ? AND data_source = ?"
data_check_result = DuckDB.execute(conn, data_check_query, [species_id, "GRI-MECH"])
data_check_df = DataFrame(data_check_result)

if size(data_check_df, 1) > 0
    println("GRI-MECH data already exists for N2")
    
    # Update reliability score
    update_data_query = """
    UPDATE thermodynamic_data 
    SET reliability_score = ?
    WHERE species_id = ? AND data_source = ?
    """
    
    DuckDB.execute(conn, update_data_query, [3.5, species_id, "GRI-MECH"])
    println("Updated reliability score for existing data")
else
    # Add to thermodynamic_data table with explicit ID
    data_query = """
    INSERT INTO thermodynamic_data (
        id, species_id, data_source, polynomial_type, 
        temperature_min, temperature_max, data_json, reliability_score
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """
    
    DuckDB.execute(conn, data_query, [
        1,  # ID
        species_id, "GRI-MECH", "nasa7", 
        200.0, 6000.0, data_json, 3.5  # Add reliability score
    ])
    
    println("Added NASA-7 polynomial data for N2")
end

# Add a data source entry if it doesn't exist
source_check_query = "SELECT name FROM data_sources WHERE name = ?"
source_check_result = DuckDB.execute(conn, source_check_query, ["GRI-MECH"])
source_check_df = DataFrame(source_check_result)

if size(source_check_df, 1) > 0
    println("GRI-MECH data source already exists")
else
    # Add to data_sources table with explicit ID
    source_query = """
    INSERT INTO data_sources (
        id, name, version, url, description, priority, reliability_score
    )
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """
    
    DuckDB.execute(conn, source_query, [
        1,  # ID
        "GRI-MECH", 
        "3.0", 
        "https://github.com/berkeleylab/gri-mech/raw/main/version30/files30/thermo30.dat",
        "GRI-Mech 3.0 Thermodynamic Database",
        1,  # Priority
        3.5  # Reliability score
    ])
    
    println("Added GRI-MECH data source")
end

# Vacuum and close the database
println("Finalizing database...")
JThermodynamicsData.vacuum_database(conn)
JThermodynamicsData.close_database(conn)

println("Minimal database initialization complete!")
println("You can now run the debug_single_species.jl script")