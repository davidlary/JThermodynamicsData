#!/usr/bin/env julia

"""
Script to update the ionic species data from various sources.
This ensures no hardcoding of thermodynamic data in the main code.
"""

using YAML
using JSON

# Determine the project root directory
const PROJECT_ROOT = dirname(dirname(abspath(@__FILE__)))

"""
    read_burcat_data(file_path::String)

Read thermodynamic data from the Burcat database file.
"""
function read_burcat_data(file_path::String)
    ionic_data = Dict()
    
    # Check if file exists
    if !isfile(file_path)
        @warn "Burcat data file not found: $file_path"
        return ionic_data
    end
    
    # Read the file and parse data for ions
    open(file_path, "r") do f
        lines = readlines(f)
        current_species = ""
        reading_coeffs = false
        coeffs_low = []
        coeffs_high = []
        
        for line in lines
            if startswith(line, "!") || isempty(strip(line))
                continue
            end
            
            if !reading_coeffs && length(line) >= 18
                # This might be a species header line
                formula = strip(line[1:18])
                if endswith(formula, "+") || endswith(formula, "-")
                    # This is an ionic species
                    current_species = strip(formula)
                    reading_coeffs = true
                    coeffs_low = []
                    coeffs_high = []
                end
            elseif reading_coeffs && length(coeffs_low) < 7
                # Reading low temperature coefficients
                if length(line) >= 75
                    for i in 1:5
                        start_pos = 1 + (i-1)*15
                        coeff_str = strip(line[start_pos:start_pos+14])
                        push!(coeffs_low, parse(Float64, coeff_str))
                    end
                end
            elseif reading_coeffs && length(coeffs_low) == 7 && length(coeffs_high) < 7
                # Reading high temperature coefficients
                if length(line) >= 75
                    for i in 1:5
                        start_pos = 1 + (i-1)*15
                        coeff_str = strip(line[start_pos:start_pos+14])
                        push!(coeffs_high, parse(Float64, coeff_str))
                    end
                end
            elseif reading_coeffs && length(coeffs_low) == 7 && length(coeffs_high) == 7
                # Finished reading this species
                ionic_data[current_species] = Dict(
                    "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
                    "coefficients" => [coeffs_low, coeffs_high],
                    "source" => "Burcat",
                    "priority" => 7,
                    "reliability_score" => 4.7
                )
                reading_coeffs = false
                current_species = ""
            end
        end
    end
    
    return ionic_data
end

"""
    read_nasa_cea_data(file_path::String)

Read thermodynamic data from the NASA CEA database file.
"""
function read_nasa_cea_data(file_path::String)
    ionic_data = Dict()
    
    # Check if file exists
    if !isfile(file_path)
        @warn "NASA CEA data file not found: $file_path"
        return ionic_data
    end
    
    # Read the file and parse data for ions (similar to Burcat but with NASA CEA format)
    open(file_path, "r") do f
        lines = readlines(f)
        current_species = ""
        reading_coeffs = false
        coeffs_low = []
        coeffs_high = []
        
        for line in lines
            if startswith(line, "#") || isempty(strip(line))
                continue
            end
            
            # NASA CEA format differs - adjust parsing as needed
            # This is a simplified example; actual parsing would be more complex
            if !reading_coeffs && length(line) >= 10
                formula = strip(line[1:10])
                if endswith(formula, "+") || endswith(formula, "-")
                    current_species = strip(formula)
                    reading_coeffs = true
                    coeffs_low = []
                    coeffs_high = []
                end
            elseif reading_coeffs
                # Read coefficient lines (simplified)
                # ...
                
                # When done reading:
                ionic_data[current_species] = Dict(
                    "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
                    "coefficients" => [coeffs_low, coeffs_high],
                    "source" => "NASA-CEA",
                    "priority" => 3,
                    "reliability_score" => 4.0
                )
                reading_coeffs = false
                current_species = ""
            end
        end
    end
    
    return ionic_data
end

"""
    merge_with_defaults(species_data::Dict)

Merge species data with defaults for missing species.
"""
function merge_with_defaults(species_data::Dict)
    # Default ionic species list if none found in data sources
    default_ionic_species = [
        "He+", "Xe+", "H+", "O+", "N+", "C+", "e-", 
        "Ar+", "Cl+", "F+", "Ne+", "Li+", "Na+", "K+", "Fe+"
    ]
    
    # Default coefficients for common ions if not found in data sources
    default_coeffs = Dict(
        "He+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 285256.641, 3.35878107],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 285256.641, 3.35878107]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "Xe+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 111732.738, 7.33529286],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 111732.738, 7.33529286]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "H+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 1529983.045, -1.14064664],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 1529983.045, -1.14064664]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "O+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 187935.285, 4.39161085],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 187935.285, 4.39161085]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "N+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 223538.655, 4.20556943],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 223538.655, 4.20556943]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "C+" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, 216142.608, 4.95832945],
                [2.500000, 0.0, 0.0, 0.0, 0.0, 216142.608, 4.95832945]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        ),
        "e-" => Dict(
            "temp_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [2.500000, 0.0, 0.0, 0.0, 0.0, -745.375, -1.17246809],
                [2.500000, 0.0, 0.0, 0.0, 0.0, -745.375, -1.17246809]
            ],
            "source" => "Default",
            "priority" => 7,
            "reliability_score" => 4.7
        )
    )
    
    # Ensure ionic_species exists
    if !haskey(species_data, "ionic_species")
        species_data["ionic_species"] = default_ionic_species
    end
    
    # Ensure nasa7_coefficients exists
    if !haskey(species_data, "nasa7_coefficients")
        species_data["nasa7_coefficients"] = Dict()
    end
    
    # Add default coefficients for missing ionic species
    for ion in species_data["ionic_species"]
        if !haskey(species_data["nasa7_coefficients"], ion) && haskey(default_coeffs, ion)
            species_data["nasa7_coefficients"][ion] = default_coeffs[ion]
        end
    end
    
    return species_data
end

"""
    update_ionic_data()

Main function to update the ionic species data file.
"""
function update_ionic_data()
    # File paths
    config_dir = joinpath(PROJECT_ROOT, "config", "data")
    output_file = joinpath(config_dir, "species_data.yaml")
    burcat_file = joinpath(PROJECT_ROOT, "data", "cache", "burcat", "burcat", "BURCAT.THR")
    nasa_cea_file = joinpath(PROJECT_ROOT, "data", "cache", "nasa-cea", "thermo.inp")
    
    # Create config directory if it doesn't exist
    mkpath(config_dir)
    
    # Initialize species data
    species_data = Dict()
    
    # Read existing file if it exists
    if isfile(output_file)
        try
            species_data = YAML.load_file(output_file)
        catch e
            @warn "Error reading existing species data file: $e"
            species_data = Dict()
        end
    end
    
    # Read data from sources
    burcat_data = read_burcat_data(burcat_file)
    nasa_cea_data = read_nasa_cea_data(nasa_cea_file)
    
    # Update nasa7_coefficients with data from sources
    if !haskey(species_data, "nasa7_coefficients")
        species_data["nasa7_coefficients"] = Dict()
    end
    
    # Add Burcat data (higher priority)
    for (species, data) in burcat_data
        species_data["nasa7_coefficients"][species] = data
    end
    
    # Add NASA CEA data (lower priority) only for species not in Burcat
    for (species, data) in nasa_cea_data
        if !haskey(species_data["nasa7_coefficients"], species)
            species_data["nasa7_coefficients"][species] = data
        end
    end
    
    # Ensure all required data is present
    species_data = merge_with_defaults(species_data)
    
    # Write updated data to file
    open(output_file, "w") do io
        YAML.write(io, species_data)
    end
    
    println("Ionic species data updated successfully at: $output_file")
    println("Found $(length(keys(species_data["nasa7_coefficients"]))) ionic/special species")
end

# Execute the update function when this script is run
if abspath(PROGRAM_FILE) == @__FILE__
    update_ionic_data()
end