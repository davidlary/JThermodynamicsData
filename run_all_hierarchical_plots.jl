#!/usr/bin/env julia

"""
Run All Hierarchical Thermodynamic Species Plots

This script runs the hierarchical plotting for all species in the database.
It leverages the simplified plotting approach that correctly implements
the hierarchical source traversal.
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using DuckDB
using DataFrames
using Printf

function main()
    # Initialize JThermodynamicsData
    config_path = joinpath(@__DIR__, "config", "settings.yaml")
    config, db = JThermodynamicsData.initialize(config_path)
    
    if db === nothing
        error("Failed to connect to database")
    end
    
    # Get all species from the database
    query = """
    SELECT name 
    FROM species 
    ORDER BY name
    """
    
    result = DuckDB.execute(db, query)
    species_df = DataFrame(result)
    
    if nrow(species_df) == 0
        error("No species found in database")
    end
    
    # Create output directory
    plots_dir = joinpath(@__DIR__, "plots_hierarchical")
    
    # Include the simple hierarchical plot module
    include(joinpath(@__DIR__, "simple_hierarchical_plot.jl"))
    
    # Count for progress reporting
    total_species = nrow(species_df)
    processed = 0
    successful = 0
    
    # Process each species
    for species in species_df.name
        processed += 1
        
        try
            println("\nProcessing $species... ($processed/$total_species)")
            plot_species_hierarchical(species, db, config, save_to_dir=plots_dir)
            successful += 1
        catch e
            println("Error processing $species: $e")
        end
        
        # Print progress
        if processed % 10 == 0
            println("Progress: $processed/$total_species species processed ($successful successful)")
        end
    end
    
    println("\nPlotting complete!")
    println("Successfully processed $successful out of $total_species species")
    println("Plots saved to: $(abspath(plots_dir))")
end

# Run the main function
main()