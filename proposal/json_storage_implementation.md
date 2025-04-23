# JSON Storage Implementation for JThermodynamicsData

## Overview

This document outlines a plan for migrating from DuckDB database storage to a JSON-per-species file-based approach for the JThermodynamicsData package.

## Storage Structure

### Directory Layout
```
data/
  species/
    H2O.json
    CO2.json
    CH4.json
    ...
  cache/
    ...  # cache files (unchanged)
  sources.json  # Source metadata and priorities
```

### JSON File Structure
Each species file (e.g., `H2O.json`) follows the schema defined in `json_storage_schema.json` and includes:

- Species metadata (name, formula, CAS, etc.)
- Data from multiple sources in a structured format
- Source priorities for hierarchical selection

## Implementation Plan

### 1. Create Core JSON I/O Functions

```julia
# src/JThermodynamicsData/database/json_storage.jl

"""
    save_species_json(species_data::Dict, output_dir::String)

Save thermodynamic data for a species to a JSON file.
"""
function save_species_json(species_data::Dict, output_dir::String)
    species_name = species_data["species_name"]
    safe_name = replace(species_name, r"[^a-zA-Z0-9_-]" => "_")
    
    # Create directory if it doesn't exist
    isdir(output_dir) || mkpath(output_dir)
    
    # Write JSON file
    open(joinpath(output_dir, "$(safe_name).json"), "w") do io
        JSON.print(io, species_data, 4)  # Pretty print with 4-space indent
    end
    
    return joinpath(output_dir, "$(safe_name).json")
end

"""
    load_species_json(species_name::String, data_dir::String)

Load thermodynamic data for a species from a JSON file.
"""
function load_species_json(species_name::String, data_dir::String)
    safe_name = replace(species_name, r"[^a-zA-Z0-9_-]" => "_")
    file_path = joinpath(data_dir, "$(safe_name).json")
    
    if !isfile(file_path)
        return nothing
    end
    
    return JSON.parsefile(file_path)
end

"""
    get_all_species(data_dir::String)

Get a list of all available species in the data directory.
"""
function get_all_species(data_dir::String)
    species_list = String[]
    
    # List all JSON files in the directory
    for file in readdir(data_dir)
        if endswith(lowercase(file), ".json")
            push!(species_list, replace(file, r"\.json$" => ""))
        end
    end
    
    return species_list
end

"""
    save_sources_metadata(sources::Dict, data_dir::String)

Save metadata about data sources to sources.json.
"""
function save_sources_metadata(sources::Dict, data_dir::String)
    # Ensure directory exists
    isdir(data_dir) || mkpath(data_dir)
    
    # Write sources file
    open(joinpath(data_dir, "sources.json"), "w") do io
        JSON.print(io, sources, 4)
    end
    
    return joinpath(data_dir, "sources.json")
end

"""
    load_sources_metadata(data_dir::String)

Load metadata about data sources from sources.json.
"""
function load_sources_metadata(data_dir::String)
    file_path = joinpath(data_dir, "sources.json")
    
    if !isfile(file_path)
        return Dict{String, Any}("sources" => Dict())
    end
    
    return JSON.parsefile(file_path)
end
```

### 2. Create Query Functions that Replace DuckDB Queries

```julia
# src/JThermodynamicsData/database/json_queries.jl

"""
    get_thermodynamic_data(species_name::String; 
                         temp_min::Float64=0.0, temp_max::Float64=1e6,
                         data_source::Union{String, Nothing}=nothing)

Get thermodynamic data for a species within a temperature range.
"""
function get_thermodynamic_data(species_name::String; 
                              temp_min::Float64=0.0, temp_max::Float64=1e6,
                              data_source::Union{String, Nothing}=nothing)
    # Load species data
    data_dir = joinpath(get_data_dir(), "species")
    species_data = load_species_json(species_name, data_dir)
    
    if species_data === nothing
        return DataFrame()
    end
    
    # Load source priorities
    sources_meta = load_sources_metadata(get_data_dir())
    source_priorities = get(sources_meta, "source_priorities", Dict())
    
    # Filter and rank sources
    filtered_sources = []
    
    for (source_id, source_data) in species_data["sources"]
        # Skip if source doesn't match requested source
        if data_source !== nothing && source_data["data_source"] != data_source
            continue
        end
        
        # Skip if temperature range doesn't overlap
        if source_data["temperature_max"] < temp_min || source_data["temperature_min"] > temp_max
            continue
        end
        
        # Get priority for this source
        priority = get(source_priorities, source_data["data_source"], 0)
        
        # Add to results
        push!(filtered_sources, (
            source_id = source_id,
            source_data = source_data,
            priority = priority,
            reliability = get(source_data, "reliability_score", 0.0)
        ))
    end
    
    # Sort by priority and reliability
    sort!(filtered_sources, by = s -> (s.priority, s.reliability), rev=true)
    
    # Convert to DataFrame
    if isempty(filtered_sources)
        return DataFrame()
    end
    
    # Create DataFrame with results
    df = DataFrame(
        species_name = species_data["species_name"],
        formula = species_data["formula"],
        cas_number = get(species_data, "cas_number", ""),
        molecular_weight = get(species_data, "molecular_weight", 0.0),
        data_source = String[],
        polynomial_type = String[],
        temperature_min = Float64[],
        temperature_max = Float64[],
        reliability_score = Float64[],
        data_json = String[],
        uncertainty_json = String[],
        priority = Int[],
        date_modified = String[]
    )
    
    # Add each source to DataFrame
    for src in filtered_sources
        source_data = src.source_data
        
        # Convert data and uncertainty to JSON strings
        data_json = JSON.json(source_data["data"])
        uncertainty_json = JSON.json(get(source_data, "uncertainty", Dict()))
        
        # Add row to DataFrame
        push!(df, (
            species_data["species_name"],
            species_data["formula"],
            get(species_data, "cas_number", ""),
            get(species_data, "molecular_weight", 0.0),
            source_data["data_source"],
            source_data["polynomial_type"],
            source_data["temperature_min"],
            source_data["temperature_max"],
            get(source_data, "reliability_score", 0.0),
            data_json,
            uncertainty_json,
            src.priority,
            get(source_data, "date_modified", "")
        ))
    end
    
    return df
end

"""
    add_species(name::String, formula::String, cas_number::String="", 
                molecular_weight::Float64=0.0, metadata::Dict=Dict())

Add a new species or update an existing one.
"""
function add_species(name::String, formula::String, cas_number::String="", 
                     molecular_weight::Float64=0.0, metadata::Dict=Dict())
    data_dir = joinpath(get_data_dir(), "species")
    isdir(data_dir) || mkpath(data_dir)
    
    # Check if species already exists
    species_data = load_species_json(name, data_dir)
    
    if species_data === nothing
        # Create new species data
        species_data = Dict{String, Any}(
            "species_name" => name,
            "formula" => formula,
            "cas_number" => cas_number,
            "molecular_weight" => molecular_weight,
            "metadata" => metadata,
            "sources" => Dict{String, Any}()
        )
    else
        # Update existing species data
        species_data["formula"] = formula
        species_data["cas_number"] = cas_number
        species_data["molecular_weight"] = molecular_weight
        species_data["metadata"] = merge(get(species_data, "metadata", Dict()), metadata)
    end
    
    # Save species data
    save_species_json(species_data, data_dir)
    
    return name
end

"""
    add_thermodynamic_data(species_name::String, data_source::String,
                          polynomial_type::String, temp_min::Float64, temp_max::Float64,
                          data::Dict, uncertainty::Dict=Dict(), reliability_score::Float64=0.0)

Add thermodynamic data for a species from a specific source.
"""
function add_thermodynamic_data(species_name::String, data_source::String,
                               polynomial_type::String, temp_min::Float64, temp_max::Float64,
                               data::Dict, uncertainty::Dict=Dict(), reliability_score::Float64=0.0)
    data_dir = joinpath(get_data_dir(), "species")
    
    # Load species data
    species_data = load_species_json(species_name, data_dir)
    
    if species_data === nothing
        @warn "Species not found: $species_name"
        return false
    end
    
    # Create source ID
    source_id = data_source
    
    # Update or add source data
    sources = get!(species_data, "sources", Dict{String, Any}())
    
    sources[source_id] = Dict{String, Any}(
        "data_source" => data_source,
        "polynomial_type" => polynomial_type,
        "temperature_min" => temp_min,
        "temperature_max" => temp_max,
        "reliability_score" => reliability_score,
        "data" => data,
        "uncertainty" => uncertainty,
        "date_added" => get(get(sources, source_id, Dict()), "date_added", string(now())),
        "date_modified" => string(now())
    )
    
    # Save updated species data
    save_species_json(species_data, data_dir)
    
    return true
end
```

### 3. Update Calculator to Use JSON Storage

```julia
# Modified version of load_polynomial_data in calculator.jl

"""
    load_polynomial_data(species_name, temperature_range, options=Dict())

Load polynomial data for thermodynamic calculations using JSON files.
"""
function load_polynomial_data(species_name, temperature_range, options=Dict())
    # Set default options
    method = get(options, "method", "hierarchical")
    source = get(options, "source", nothing)
    enable_uncertainty = get(options, "uncertainty", true)
    
    # Query the data from JSON storage
    if source !== nothing
        df = get_thermodynamic_data(species_name; 
            temp_min=temperature_range[1], 
            temp_max=temperature_range[2],
            data_source=source)
    else
        df = get_thermodynamic_data(species_name; 
            temp_min=temperature_range[1], 
            temp_max=temperature_range[2])
    end
    
    if size(df, 1) == 0
        if method == "hierarchical" && get(options, "enable_fallback", true)
            # Try theoretical calculation as fallback
            return calculate_theoretical_properties(species_name, temperature_range)
        else
            error("No data found for species $(species_name) in temperature range $(temperature_range)")
        end
    end
    
    # Get first row (highest priority/reliability)
    row = df[1, :]
    
    # Parse the data JSON
    data = JSON.parse(row.data_json)
    
    # Extract coefficients and temperature ranges
    temp_ranges = [
        Float64.(data["low_temp"]["range"]),
        Float64.(data["high_temp"]["range"])
    ]
    
    coeffs = [
        Float64.(data["low_temp"]["coefficients"]),
        Float64.(data["high_temp"]["coefficients"])
    ]
    
    # Parse uncertainty information
    uncertainty_data = JSON.parse(row.uncertainty_json)
    uncertainty_value = get(uncertainty_data, "value", 0.1)
    
    # Create and return the polynomial
    return ThermodynamicPolynomial(
        row.species_name,
        row.formula,
        row.data_source,
        row.priority,
        row.reliability_score,
        row.polynomial_type,
        temp_ranges,
        coeffs,
        uncertainty_value
    )
end
```

### 4. Migration Utility

```julia
# scripts/migrate_to_json.jl

using JThermodynamicsData
using DuckDB
using DataFrames
using JSON
using Dates

"""
    migrate_database_to_json(db_path::String, output_dir::String)

Migrate data from DuckDB database to JSON files.
"""
function migrate_database_to_json(db_path::String, output_dir::String)
    # Connect to database
    conn = DBInterface.connect(DuckDB.DB, db_path)
    
    # Create output directories
    species_dir = joinpath(output_dir, "species")
    mkpath(species_dir)
    
    # Migrate data sources
    migrate_data_sources(conn, output_dir)
    
    # Get all species
    species_result = DuckDB.execute(conn, "SELECT id, name, formula, cas_number, molecular_weight, metadata_json FROM species")
    species_df = DataFrame(species_result)
    
    # Process each species
    for row in eachrow(species_df)
        species_id = row.id
        species_name = row.name
        
        println("Migrating species: $species_name")
        
        # Create species data structure
        species_data = Dict{String, Any}(
            "species_name" => species_name,
            "formula" => row.formula,
            "cas_number" => row.cas_number,
            "molecular_weight" => row.molecular_weight,
            "sources" => Dict{String, Any}()
        )
        
        # Add metadata if available
        if !ismissing(row.metadata_json) && row.metadata_json != "{}"
            try
                species_data["metadata"] = JSON.parse(row.metadata_json)
            catch e
                @warn "Error parsing metadata for $species_name: $e"
                species_data["metadata"] = Dict()
            end
        end
        
        # Get thermodynamic data for this species
        thermo_query = """
        SELECT 
            data_source,
            polynomial_type,
            temperature_min,
            temperature_max,
            reliability_score,
            data_json,
            uncertainty_json,
            date_added,
            date_modified
        FROM 
            thermodynamic_data
        WHERE 
            species_id = ?
        """
        
        thermo_result = DuckDB.execute(conn, thermo_query, [species_id])
        thermo_df = DataFrame(thermo_result)
        
        # Process each data source
        for trow in eachrow(thermo_df)
            source_id = trow.data_source
            
            # Parse data and uncertainty JSON
            data = JSON.parse(trow.data_json)
            uncertainty = Dict()
            if !ismissing(trow.uncertainty_json) && trow.uncertainty_json != "{}"
                try
                    uncertainty = JSON.parse(trow.uncertainty_json)
                catch e
                    @warn "Error parsing uncertainty for $species_name/$source_id: $e"
                end
            end
            
            # Add source data
            species_data["sources"][source_id] = Dict{String, Any}(
                "data_source" => trow.data_source,
                "polynomial_type" => trow.polynomial_type,
                "temperature_min" => trow.temperature_min,
                "temperature_max" => trow.temperature_max,
                "reliability_score" => trow.reliability_score,
                "data" => data,
                "uncertainty" => uncertainty,
                "date_added" => string(trow.date_added),
                "date_modified" => string(trow.date_modified)
            )
        end
        
        # Save species data to JSON file
        save_species_json(species_data, species_dir)
    end
    
    println("Migration complete!")
end

"""
    migrate_data_sources(conn::DuckDB.DB, output_dir::String)

Migrate data sources metadata to JSON.
"""
function migrate_data_sources(conn::DuckDB.DB, output_dir::String)
    # Get data sources
    sources_result = DuckDB.execute(conn, """
    SELECT 
        name, 
        version, 
        url, 
        description,
        priority,
        last_updated,
        metadata_json
    FROM 
        data_sources
    ORDER BY 
        priority DESC
    """)
    
    sources_df = DataFrame(sources_result)
    
    # Create sources metadata structure
    sources_meta = Dict{String, Any}(
        "sources" => Dict{String, Any}(),
        "source_priorities" => Dict{String, Int}()
    )
    
    # Process each source
    for row in eachrow(sources_df)
        source_name = row.name
        
        # Add source metadata
        sources_meta["sources"][source_name] = Dict{String, Any}(
            "name" => source_name,
            "version" => row.version,
            "url" => row.url,
            "description" => row.description,
            "last_updated" => string(row.last_updated)
        )
        
        # Add source priority
        sources_meta["source_priorities"][source_name] = row.priority
        
        # Add additional metadata if available
        if !ismissing(row.metadata_json) && row.metadata_json != "{}"
            try
                metadata = JSON.parse(row.metadata_json)
                sources_meta["sources"][source_name]["metadata"] = metadata
            catch e
                @warn "Error parsing metadata for source $source_name: $e"
            end
        end
    end
    
    # Save sources metadata
    open(joinpath(output_dir, "sources.json"), "w") do io
        JSON.print(io, sources_meta, 4)
    end
    
    println("Migrated $(length(sources_meta["sources"])) data sources")
end

# Run migration if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        println("Usage: julia migrate_to_json.jl <database_path> <output_directory>")
        exit(1)
    end
    
    db_path = ARGS[1]
    output_dir = ARGS[2]
    
    migrate_database_to_json(db_path, output_dir)
end
```

## Integration with Existing Code

To integrate with the existing codebase:

1. Implement the JSON storage and query modules as outlined above
2. Create a compatibility layer that makes the JSON storage look like the DuckDB interface:
   - Implement database connection as a path to the species directory
   - Replace existing query functions with JSON-based versions
   - Update the calculator to work with the new JSON backend

## Pros and Cons Analysis

### Pros of JSON-per-species Approach

1. **Simplicity**: Removes database dependency, making the package easier to install and use
2. **Incremental Updates**: Files can be updated independently as new data becomes available
3. **Version Control Friendly**: Individual JSON files are easier to track in version control
4. **Direct Access**: Data can be accessed directly without SQL queries
5. **Portability**: JSON files can be easily moved between systems
6. **Read/Write Performance**: Better for read-heavy workloads (thermodynamic lookups)
7. **Transparency**: Data structure is visible and human-readable
8. **Hierarchical Selection**: Can still be implemented effectively with source priorities

### Cons of JSON-per-species Approach

1. **No SQL Capabilities**: Complex filtering and joins must be implemented in Julia code
2. **No Transactional Safety**: No ACID guarantees for concurrent writes
3. **No Indexes**: Must implement search capabilities manually
4. **Memory Usage**: May require loading more data into memory
5. **Schema Enforcement**: Requires manual validation of JSON schema
6. **Cross-Species Queries**: More complex and potentially slower
7. **Bulk Operations**: Updating many species at once becomes more complex

### Implementation Difficulty

The transition difficulty is **medium**:

1. **Core Functionality**: Most code changes are isolated to database interactions
2. **Migration Utility**: A one-time migration script is needed
3. **API Compatibility**: The external API can remain mostly unchanged
4. **Testing**: Thorough testing will be required to ensure equivalent behavior

## Conclusion

The JSON-per-species approach offers a good balance of simplicity, flexibility, and performance for the JThermodynamicsData package. While it sacrifices some features of a database, the benefits of simpler installation and more transparent data storage make it an attractive option. The implementation effort is reasonable and can be done incrementally.

The proposed implementation maintains compatibility with the existing API while moving to a more flexible storage backend that should better serve the package's long-term needs.