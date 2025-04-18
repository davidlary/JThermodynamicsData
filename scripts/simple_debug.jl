#!/usr/bin/env julia

# Very simple script to calculate thermodynamic properties for N2 at a specific temperature
# This is a simplified version of debug_single_species.jl

using Pkg
Pkg.activate(".")

println("Loading JThermodynamicsData module...")
using JThermodynamicsData
using YAML
using DuckDB
using DataFrames
using JSON

println("Initializing configuration and database connection...")
config_path = joinpath(dirname(@__DIR__), "config", "settings.yaml")
config = JThermodynamicsData.load_config(config_path)
conn = JThermodynamicsData.init_database(config)

# Set up a single temperature to test
species_name = "N2"
temperature = 298.15  # Standard temperature in K

println("Processing thermodynamic properties for $species_name at $temperature K...")

try
    # Check if species exists in the database
    species_id = JThermodynamicsData.get_species_id(conn, species_name)
    
    if species_id === nothing
        println("Species not found in database")
    else
        println("Species found in database with ID: $species_id")
        
        # Test direct query to database for this species
        query_result = DuckDB.execute(conn, "SELECT * FROM species WHERE name = ?", [species_name])
        species_data = DataFrame(query_result)
        println("Species data: ", species_data)
        
        # Get thermodynamic data sources for this species
        data_query = "SELECT * FROM thermodynamic_data WHERE species_id = ?"
        data_result = DuckDB.execute(conn, data_query, [species_id])
        data_df = DataFrame(data_result)
        println("\nThermodynamic data available:")
        println(data_df)
        
        # Update species with necessary metadata
        println("\nUpdating species metadata...")
        metadata = Dict(
            "molecular_data" => Dict(
                "molecular_weight" => 28.0134,
                "symmetry_number" => 2,
                "spin_multiplicity" => 1,
                "rotational_constants" => [1.998],  # cm^-1
                "vibrational_frequencies" => [2358.0]  # cm^-1
            )
        )
        metadata_json = JSON.json(metadata)
        
        update_query = "UPDATE species SET metadata_json = ? WHERE id = ?"
        DuckDB.execute(conn, update_query, [metadata_json, species_id])
        println("Updated metadata for N2")
        
        # Update data sources to include priority information
        priority_query = """
        SELECT name FROM data_sources WHERE name = 'GRI-MECH'
        """
        source_result = DuckDB.execute(conn, priority_query)
        source_df = DataFrame(source_result)
        
        if size(source_df, 1) == 0
            source_query = """
            INSERT INTO data_sources (id, name, version, url, priority, reliability_score) 
            VALUES (1, 'GRI-MECH', '3.0', 'https://github.com/berkeleylab/gri-mech/raw/main/version30/files30/thermo30.dat', 1, 3.5)
            """
            DuckDB.execute(conn, source_query)
            println("Added GRI-MECH source with priority data")
        end
        
        # Fix issue with constants by monkeypatching
        println("\nApplying workaround for constants...")
        JThermodynamicsData.h = JThermodynamicsData.H  # Use uppercase constant
        JThermodynamicsData.c = JThermodynamicsData.C
        JThermodynamicsData.kb = JThermodynamicsData.KB
        JThermodynamicsData.na = JThermodynamicsData.NA
        println("Applied lowercase constant workarounds")
    end
    
    # Try to get properties directly from the database
    println("\nAttempting to query properties directly...")
    
    # First, get data from database for the specified temperature
    @info "Getting thermodynamic data for $species_name at $temperature K"
    
    # Use progressively_refine_thermodynamic_data which is the core function
    result, sources = JThermodynamicsData.progressively_refine_thermodynamic_data(
        conn, species_name, temperature, config
    )
    
    println("\nResults using progressively_refine_thermodynamic_data:")
    println("Data source: $(result["data_source"])")
    println("Properties:")
    for (prop, data) in result["properties"]
        println("  $prop: $(data["value"]) ± $(data["uncertainty"]) $(data["units"])")
    end
    
    println("\nNumber of data sources used: $(length(sources))")
    for (i, source) in enumerate(sources)
        println("Source $i: $(source["data_source"])")
    end
    
catch e
    println("ERROR: $e")
    println(stacktrace(catch_backtrace()))
finally
    # Close the database connection
    println("\nClosing database connection...")
    JThermodynamicsData.close_database(conn)
end

println("Debug run complete!")