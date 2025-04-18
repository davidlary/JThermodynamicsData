"""
Functions for caching thermodynamic data and calculations.
"""

"""
    get_cache_key(species_name::String, temperature::Float64, property_type::String)

Generate a cache key for a specific species, temperature, and property.
"""
function get_cache_key(species_name::String, temperature::Float64, property_type::String)
    return "$(lowercase(species_name))_$(round(temperature, digits=2))_$(property_type)"
end

"""
    store_in_cache(conn::DuckDB.DB, species_id::Int, temperature_min::Float64, 
                 temperature_max::Float64, cache_type::String, data::Dict)

Store data in the cache table.
"""
function store_in_cache(conn::DuckDB.DB, species_id::Int, temperature_min::Float64, 
                      temperature_max::Float64, cache_type::String, data::Dict)
    # Convert data to JSON
    data_json = JSON.json(data)
    
    # Check if cache entry already exists
    result = DuckDB.execute(conn, """
    SELECT id FROM data_cache
    WHERE species_id = ? AND temperature_min = ? AND temperature_max = ? AND cache_type = ?
    LIMIT 1
    """, [species_id, temperature_min, temperature_max, cache_type])
    
    df = DataFrame(result)
    
    if size(df, 1) > 0
        # Update existing cache entry
        record_id = df[1, :id]
        
        stmt = DuckDB.prepare(conn, """
        UPDATE data_cache
        SET cache_data = ?, created_at = CURRENT_TIMESTAMP
        WHERE id = ?
        """)
        
        DuckDB.execute(stmt, [data_json, record_id])
        return record_id
    else
        # Insert new cache entry
        stmt = DuckDB.prepare(conn, """
        INSERT INTO data_cache (species_id, temperature_min, temperature_max, cache_type, cache_data)
        VALUES (?, ?, ?, ?, ?)
        RETURNING id
        """)
        
        result = DuckDB.execute(stmt, [species_id, temperature_min, temperature_max, cache_type, data_json])
        df = DataFrame(result)
        return df[1, :id]
    end
end

"""
    get_from_cache(conn::DuckDB.DB, species_id::Int, temperature::Real, cache_type::String)

Retrieve data from the cache table.
"""
function get_from_cache(conn::DuckDB.DB, species_id::Int, temperature::Real, cache_type::String)
    # Convert temperature to Float64 to ensure type stability
    temp = Float64(temperature)
    
    # Query cache table for entries that cover the requested temperature
    result = DuckDB.execute(conn, """
    SELECT id, cache_data
    FROM data_cache
    WHERE species_id = ? AND temperature_min <= ? AND temperature_max >= ? AND cache_type = ?
    ORDER BY created_at DESC
    LIMIT 1
    """, [species_id, temp, temp, cache_type])
    
    df = DataFrame(result)
    
    if size(df, 1) > 0
        # Parse cached data
        try
            cache_data = JSON.parse(df[1, :cache_data])
            return cache_data
        catch e
            @warn "Failed to parse cached data: $e"
            return nothing
        end
    else
        return nothing
    end
end

"""
    clear_cache(conn::DuckDB.DB, species_id::Union{Int, Nothing}=nothing, 
              cache_type::Union{String, Nothing}=nothing)

Clear cache entries for a species and/or cache type.
"""
function clear_cache(conn::DuckDB.DB, species_id::Union{Int, Nothing}=nothing, 
                   cache_type::Union{String, Nothing}=nothing)
    if species_id === nothing && cache_type === nothing
        # Clear all cache
        DuckDB.execute(conn, "DELETE FROM data_cache")
        @info "Cleared all cache entries"
    elseif species_id !== nothing && cache_type === nothing
        # Clear cache for species
        stmt = DuckDB.prepare(conn, "DELETE FROM data_cache WHERE species_id = ?")
        DuckDB.execute(stmt, [species_id])
        @info "Cleared cache entries for species ID: $species_id"
    elseif species_id === nothing && cache_type !== nothing
        # Clear cache by type
        stmt = DuckDB.prepare(conn, "DELETE FROM data_cache WHERE cache_type = ?")
        DuckDB.execute(stmt, [cache_type])
        @info "Cleared cache entries of type: $cache_type"
    else
        # Clear cache for species and type
        stmt = DuckDB.prepare(conn, "DELETE FROM data_cache WHERE species_id = ? AND cache_type = ?")
        DuckDB.execute(stmt, [species_id, cache_type])
        @info "Cleared cache entries for species ID: $species_id and type: $cache_type"
    end
    
    return true
end

"""
    cache_frequent_species(conn::DuckDB.DB, config::Dict)

Pre-cache data for frequently accessed species.
"""
function cache_frequent_species(conn::DuckDB.DB, config::Dict)
    # Get list of most frequently accessed species
    result = DuckDB.execute(conn, """
    SELECT species_id, COUNT(*) as query_count
    FROM data_cache
    GROUP BY species_id
    ORDER BY query_count DESC
    LIMIT 50
    """)
    
    df = DataFrame(result)
    
    if size(df, 1) == 0
        @info "No frequent species to cache"
        return
    end
    
    # Default temperature range
    temp_range = config["general"]["temperature_range"]
    temp_min = temp_range[1]
    temp_max = temp_range[2]
    
    # Precalculate temperatures at regular intervals
    temp_step = 100.0  # 100K intervals
    
    precalc_temps = temp_min:temp_step:temp_max
    
    # For each species, precalculate properties at standard temperatures
    for row in eachrow(df)
        species_id = row.species_id
        
        # Get species name
        species_result = DuckDB.execute(conn, "SELECT name FROM species WHERE id = ?", [species_id])
        species_df = DataFrame(species_result)
        
        if size(species_df, 1) == 0
            continue
        end
        
        species_name = species_df[1, :name]
        
        @info "Pre-caching data for frequent species: $species_name"
        
        # Query properties at each temperature
        for temp in precalc_temps
            # Skip if already in cache
            cached = get_from_cache(conn, species_id, temp, "properties")
            
            if cached !== nothing
                continue
            end
            
            # Query properties
            try
                data = query_properties(conn, species_name, temp, config)
                
                # Store in cache
                store_in_cache(conn, species_id, temp, temp, "properties", data)
            catch e
                @warn "Failed to cache properties for $species_name at $temp K: $e"
            end
        end
    end
    
    @info "Completed pre-caching for frequent species"
    return true
end