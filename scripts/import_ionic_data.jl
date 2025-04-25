#!/usr/bin/env julia

"""
Import Ionic Species Data to Database

This script imports the thermodynamic data for ionic species from species_data.yaml
directly into the database. This ensures all data is stored in a single location
(the database) rather than being split between files and database.

After running this script, the database will be the single source of truth for
all thermodynamic data, including ionic species.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using JSON
using DuckDB
using DataFrames
using Dates

# Function to load species_data.yaml
function load_species_data()
    project_dir = dirname(dirname(@__FILE__))
    species_data_path = joinpath(project_dir, "config", "species_data.yaml")
    
    if !isfile(species_data_path)
        error("species_data.yaml not found at $(species_data_path)")
    end
    
    try
        return YAML.load_file(species_data_path)
    catch e
        error("Error parsing species_data.yaml: $(e)")
    end
end

# Function to import thermodynamic data for a species
function import_species_thermodynamic_data(conn::DuckDB.DB, species_name::String, 
                                         data::Dict, data_source::String, priority::Int)
    # Extract basic properties
    formula = species_name  # For ionic species, the name is typically the formula
    
    # Add the species to the database if it doesn't exist
    species_id = JThermodynamicsData.add_species(conn, species_name, formula)
    
    # Extract temperature ranges and coefficients
    temp_ranges = data["temperature_ranges"]
    temp_min = temp_ranges[1][1]
    temp_max = temp_ranges[end][2]
    
    # Create data structure in NASA-7 format
    properties = data["properties"]
    uncertainty = get(data, "uncertainty", 0.05)  # Default 5% uncertainty
    reliability_score = 4.0 + (14 - priority) * 0.07  # Higher priority = higher reliability
    citation = get(data, "citation", "")
    
    # Create NASA-7 polynomial data format
    nasa7_data = Dict(
        "polynomial_type" => "nasa7",
        "temperature_ranges" => temp_ranges,
        "source" => data_source,
        "citation" => citation,
        "data" => Dict(
            "properties" => properties
        )
    )
    
    # Store uncertainty in a separate structure
    uncertainty_data = Dict(
        "method" => "literature",
        "value" => uncertainty,
        "citation" => citation
    )
    
    # Add to database
    JThermodynamicsData.add_thermodynamic_data(
        conn, species_id, data_source, "nasa7",
        temp_min, temp_max, nasa7_data, uncertainty_data, reliability_score
    )
    
    return species_id
end

# Function to import NASA-7 coefficients for a species
function import_nasa7_coefficients(conn::DuckDB.DB, species_name::String, 
                                 coeffs::Vector, data_source::String, priority::Int)
    # Extract basic properties
    formula = species_name  # For ionic species, the name is typically the formula
    
    # Add the species to the database if it doesn't exist
    species_id = JThermodynamicsData.add_species(conn, species_name, formula)
    
    # Standard temperature ranges for NASA-7
    temp_ranges = [[200.0, 1000.0], [1000.0, 6000.0]]
    temp_min = temp_ranges[1][1]
    temp_max = temp_ranges[2][2]
    
    # Set reliability based on priority
    reliability_score = 4.0 + (14 - priority) * 0.07  # Higher priority = higher reliability
    
    # Create NASA-7 polynomial data
    nasa7_data = Dict(
        "polynomial_type" => "nasa7",
        "temperature_min" => temp_min,
        "temperature_max" => temp_max,
        "data" => Dict(
            "low_temp" => Dict(
                "range" => temp_ranges[1],
                "coefficients" => coeffs[1]
            ),
            "high_temp" => Dict(
                "range" => temp_ranges[2],
                "coefficients" => coeffs[2]
            )
        )
    )
    
    # Store uncertainty (default 5%)
    uncertainty_data = Dict(
        "method" => "literature",
        "value" => 0.05
    )
    
    # Add to database
    JThermodynamicsData.add_thermodynamic_data(
        conn, species_id, data_source, "nasa7",
        temp_min, temp_max, nasa7_data, uncertainty_data, reliability_score
    )
    
    return species_id
end

# Main function to import all data
function main()
    # Set up logging
    log_dir = joinpath(dirname(dirname(@__FILE__)), "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "import_ionic_data_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start logging
    @info "Starting import of ionic species data"
    @info "=================================="
    
    # Connect to database
    db_path = joinpath(dirname(dirname(@__FILE__)), "data", "thermodynamics.duckdb")
    conn = JThermodynamicsData.get_connection(db_path)
    
    # Load species data
    @info "Loading species_data.yaml"
    species_data = load_species_data()
    
    # Import special_species data
    if haskey(species_data, "special_species")
        @info "Importing special species thermodynamic data"
        special_species_data = species_data["special_species"]
        
        for (species_name, data) in special_species_data
            @info "  Importing $(species_name)"
            
            # Use a theoretical data source with high priority
            data_source = "THEORETICAL_IONIC"
            priority = 4  # Highest theoretical priority per README
            
            # Import the data
            import_species_thermodynamic_data(conn, species_name, data, data_source, priority)
        end
        
        @info "Imported $(length(special_species_data)) special species"
    end
    
    # Import additional_species_data (NASA-7 coefficients)
    if haskey(species_data, "additional_species_data")
        @info "Importing additional species NASA-7 coefficients"
        additional_data = species_data["additional_species_data"]
        
        for (species_name, coeffs) in additional_data
            @info "  Importing NASA-7 coefficients for $(species_name)"
            
            # Use a dedicated data source for NASA-7 coefficients
            data_source = "NASA7_IONIC"
            priority = 7  # NASA-CEA priority per README
            
            # Import the coefficients
            import_nasa7_coefficients(conn, species_name, coeffs, data_source, priority)
        end
        
        @info "Imported NASA-7 coefficients for $(length(additional_data)) species"
    end
    
    # Close database connection
    DuckDB.close!(conn)
    
    @info "Ionic species data import complete!"
    @info "Now all thermodynamic data is stored in the database."
    @info "For consistency, you may want to remove species_data.yaml or rename it to species_data.yaml.backup"
end

# Run the main function
main()