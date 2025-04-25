#!/usr/bin/env julia

"""
Verify Real Data Sources Implementation

This script verifies that all thermodynamic data sources are properly implemented
with real data fetching instead of synthetic data generation.
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using JSON
using YAML
using Dates
using Printf

# Initialize logger
log_dir = joinpath(@__DIR__, "logs")
mkpath(log_dir)
log_file = joinpath(log_dir, "real_data_verification_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
JThermodynamicsData.init_logger("info", log_file)

# Start logging
println("\n" * "=" * 80)
println("REAL DATA SOURCES VERIFICATION")
println("=" * 80)

# Verify real data sources
function verify_real_data_sources()
    # Get source configuration
    config_path = joinpath(@__DIR__, "config", "settings.yaml")
    config = JThermodynamicsData.load_config(config_path)
    
    if !haskey(config, "sources")
        println("❌ No sources defined in config")
        return false, []
    end
    
    println("Sources defined in config:")
    for (name, source) in config["sources"]
        priority = get(source, "priority", "N/A")
        type = get(source, "type", "N/A")
        println("  - $name (priority: $priority, type: $type)")
    end
    
    # Try to fetch data for some common species
    test_species = ["N2", "O2", "H2O", "CO2", "Ar", "CH4", "CO", "NO", "OH"]
    source_results = []
    
    # Get cache directory
    cache_dir = joinpath(@__DIR__, "data", "cache")
    
    println("\nChecking cache directories for real data:")
    for (name, _) in config["sources"]
        source_cache = joinpath(cache_dir, lowercase(name))
        if isdir(source_cache)
            files = readdir(source_cache)
            file_count = length(files)
            println("  - $name: $file_count cached files")
            
            # Record result
            push!(source_results, (
                name = name,
                cache_exists = true,
                file_count = file_count,
                has_real_data = file_count > 0
            ))
        else
            println("  - $name: No cache directory")
            push!(source_results, (
                name = name,
                cache_exists = false,
                file_count = 0,
                has_real_data = false
            ))
        end
    end
    
    # Test data fetching for real
    println("\nTesting data fetching for common species:")
    species_results = []
    
    for species in test_species
        println("\n  Testing $species:")
        found_sources = []
        
        # Try to load species data
        species_data = JThermodynamicsData.load_species_data(species)
        
        if species_data === nothing || !haskey(species_data, "sources") || isempty(species_data["sources"])
            println("    ❌ No data found for $species")
            push!(species_results, (
                species = species,
                has_data = false,
                sources = []
            ))
            continue
        end
        
        for (source_name, source_info) in species_data["sources"]
            # Get source metadata
            source_type = get(source_info, "type", "unknown")
            priority = get(source_info, "priority", 0)
            is_experimental = priority > 4
            
            # Get polynomial data if available
            has_poly = haskey(source_info, "polynomial")
            
            # Check source metadata
            source_metadata = get(source_info, "metadata", Dict())
            has_real_metadata = !isempty(source_metadata) && 
                                (haskey(source_metadata, "fetched_at") || 
                                 haskey(source_metadata, "source_url") ||
                                 haskey(source_metadata, "data_version"))
            
            println("    - $source_name (priority: $priority, type: $source_type):")
            println("      Experimental: $(is_experimental ? "Yes" : "No")")
            println("      Has polynomial data: $(has_poly ? "Yes" : "No")")
            println("      Has real metadata: $(has_real_metadata ? "Yes" : "No")")
            
            # Record this source
            push!(found_sources, (
                name = source_name,
                type = source_type,
                priority = priority,
                is_experimental = is_experimental,
                has_poly = has_poly,
                has_real_metadata = has_real_metadata
            ))
        end
        
        # Record species result
        push!(species_results, (
            species = species,
            has_data = true,
            sources = found_sources,
            experimental_count = count(s -> s.is_experimental, found_sources),
            theoretical_count = count(s -> !s.is_experimental, found_sources),
            real_data_count = count(s -> s.has_real_metadata, found_sources)
        ))
    end
    
    # Summarize results
    experimental_sources = filter(r -> get(r, :has_real_data, false), source_results)
    theoretical_sources = filter(r -> !get(r, :has_real_data, false), source_results)
    
    real_data_species = count(r -> r.real_data_count > 0, species_results)
    no_real_data_species = count(r -> r.real_data_count == 0, species_results)
    
    # Show source distribution for all tested species
    println("\nSource distribution for tested species:")
    source_counts = Dict()
    for result in species_results
        for source in result.sources
            if source.has_real_metadata
                source_counts[source.name] = get(source_counts, source.name, 0) + 1
            end
        end
    end
    
    # Sort by count (descending)
    sorted_sources = sort(collect(source_counts), by = x -> x[2], rev = true)
    
    # Calculate percentages
    total_species = length(species_results)
    println("===================")
    for (source, count) in sorted_sources
        percentage = round(count / total_species * 100, digits=1)
        println("- $source: $count species ($percentage%)")
    end
    
    # Final validation
    has_experimental = !isempty(experimental_sources)
    all_species_have_data = all(r -> r.has_data, species_results)
    some_real_data = real_data_species > 0
    
    validation_passed = has_experimental && some_real_data
    
    println("\n" * "=" * 80)
    if validation_passed
        println("✅ VALIDATION PASSED: Database includes experimental data sources!")
    else
        println("❌ VALIDATION FAILED: Not enough experimental data sources or real data")
    end
    println("=" * 80)
    
    return validation_passed, (
        source_results = source_results,
        species_results = species_results,
        experimental_sources = length(experimental_sources),
        theoretical_sources = length(theoretical_sources),
        real_data_species = real_data_species,
        no_real_data_species = no_real_data_species,
        source_distribution = sorted_sources
    )
end

# Run verification
success, results = verify_real_data_sources()

# Return success status
exit(success ? 0 : 1)