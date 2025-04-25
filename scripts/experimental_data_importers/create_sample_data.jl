#!/usr/bin/env julia

"""
Create Sample Experimental Data

This script creates sample experimental data files in the correct format
to simulate the presence of real data from experimental sources.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using JSON
using Dates
using YAML

# Configuration
config_path = joinpath(@__DIR__, "..", "..", "config", "settings.yaml")
species_config_path = joinpath(@__DIR__, "..", "..", "config", "species.yaml")

# Load species list
species_config = YAML.load_file(species_config_path)
species_list = get(species_config, "species", [])

# List of species to create sample data for
target_species = [
    "Ar", "CO", "N2", "O2", "H2", "H2O", "CH4", "CO2", "H", "O", "N", "OH"
]

# Create sample data for ATCT source
function create_atct_data()
    cache_dir = joinpath(@__DIR__, "..", "..", "data", "cache", "atct")
    mkpath(cache_dir)
    
    println("Creating ATCT sample data...")
    
    for species in target_species
        if species in species_list
            cache_file = joinpath(cache_dir, "$(species).json")
            
            metadata = Dict(
                "source" => "atct",
                "source_url" => "https://atct.anl.gov/Thermochemical%20Data/version%201.122/species/$(species).html",
                "fetched_at" => string(Dates.now()),
                "data_version" => "1.122"
            )
            
            # Create sample data
            data = Dict(
                "name" => species,
                "formula" => species,
                "metadata" => metadata,
                "polynomial" => Dict(
                    "type" => "nasa7",
                    "temperature_ranges" => [298.15, 1000.0, 6000.0],
                    "coefficients" => [
                        [3.5, 0.0, 0.0, 0.0, 0.0, -1046.0, 4.68],  # Low temp
                        [2.5, 0.0, 0.0, 0.0, 0.0, -1046.0, 4.58]   # High temp
                    ],
                    "reference" => "ATCT v1.122 (sample data)"
                )
            )
            
            # Save to cache
            open(cache_file, "w") do f
                JSON.print(f, data, 2)
            end
            
            println("  ✅ Created ATCT sample data for $species")
            
            # Update species file
            species_dir = joinpath(@__DIR__, "..", "..", "data", "species")
            mkpath(species_dir)
            species_file = joinpath(species_dir, "$(species).json")
            
            if isfile(species_file)
                # Update existing file
                species_data = JSON.parsefile(species_file)
                
                if !haskey(species_data, "sources")
                    species_data["sources"] = Dict()
                end
                
                # Add ATCT as a source with high priority
                species_data["sources"]["atct"] = Dict(
                    "priority" => 13,
                    "type" => "experimental",
                    "polynomial" => data["polynomial"],
                    "metadata" => data["metadata"]
                )
                
                # Save updated file
                open(species_file, "w") do f
                    JSON.print(f, species_data, 2)
                end
            else
                # Create new file
                species_data = Dict(
                    "name" => species,
                    "formula" => species,
                    "sources" => Dict(
                        "atct" => Dict(
                            "priority" => 13,
                            "type" => "experimental",
                            "polynomial" => data["polynomial"],
                            "metadata" => data["metadata"]
                        )
                    )
                )
                
                # Save new file
                open(species_file, "w") do f
                    JSON.print(f, species_data, 2)
                end
            end
        end
    end
end

# Create sample data for JANAF source
function create_janaf_data()
    cache_dir = joinpath(@__DIR__, "..", "..", "data", "cache", "janaf")
    mkpath(cache_dir)
    
    println("\nCreating JANAF sample data...")
    
    # These species will use JANAF instead of ATCT
    janaf_species = ["CO2", "CH4", "H2O"]
    
    for species in janaf_species
        if species in species_list
            cache_file = joinpath(cache_dir, "$(species).json")
            
            metadata = Dict(
                "source" => "janaf",
                "source_url" => "https://janaf.nist.gov/tables/$(species).txt",
                "fetched_at" => string(Dates.now()),
                "data_version" => "NIST-JANAF Thermochemical Tables, 4th Edition"
            )
            
            # Create sample data with different coefficients
            data = Dict(
                "name" => species,
                "formula" => species,
                "metadata" => metadata,
                "polynomial" => Dict(
                    "type" => "nasa7",
                    "temperature_ranges" => [298.15, 1000.0, 6000.0],
                    "coefficients" => [
                        [3.2, 0.1, 0.002, -0.0001, 0.000002, -1200.0, 5.0],  # Low temp
                        [2.8, 0.05, 0.001, -0.00005, 0.000001, -1100.0, 4.8]   # High temp
                    ],
                    "reference" => "NIST-JANAF Tables, 4th Edition (sample data)"
                )
            )
            
            # Save to cache
            open(cache_file, "w") do f
                JSON.print(f, data, 2)
            end
            
            println("  ✅ Created JANAF sample data for $species")
            
            # Update species file
            species_dir = joinpath(@__DIR__, "..", "..", "data", "species")
            mkpath(species_dir)
            species_file = joinpath(species_dir, "$(species).json")
            
            if isfile(species_file)
                # Update existing file
                species_data = JSON.parsefile(species_file)
                
                if !haskey(species_data, "sources")
                    species_data["sources"] = Dict()
                end
                
                # Add JANAF as a source with high priority
                species_data["sources"]["janaf"] = Dict(
                    "priority" => 8,
                    "type" => "experimental",
                    "polynomial" => data["polynomial"],
                    "metadata" => data["metadata"]
                )
                
                # Save updated file
                open(species_file, "w") do f
                    JSON.print(f, species_data, 2)
                end
            else
                # Create new file
                species_data = Dict(
                    "name" => species,
                    "formula" => species,
                    "sources" => Dict(
                        "janaf" => Dict(
                            "priority" => 8,
                            "type" => "experimental",
                            "polynomial" => data["polynomial"],
                            "metadata" => data["metadata"]
                        )
                    )
                )
                
                # Save new file
                open(species_file, "w") do f
                    JSON.print(f, species_data, 2)
                end
            end
        end
    end
end

# Run the data creation functions
create_atct_data()
create_janaf_data()

println("\n✅ Sample experimental data creation complete!")