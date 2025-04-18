"""
Functions for querying thermodynamic data.
"""

"""
    query_species(conn::DuckDB.DB, species_name::String, config::Dict)

Query information about a species.
"""
function query_species(conn::DuckDB.DB, species_name::String, config::Dict)
    # Try to find species by name
    species_id = get_species_id(conn, species_name)
    
    if species_id === nothing
        # Try to search for similar species
        similar_species = search_species(conn, species_name, limit=5)
        
        if size(similar_species, 1) > 0
            suggestions = similar_species.name
            error("Species not found: $species_name. Similar species: $(join(suggestions, ", "))")
        else
            error("Species not found: $species_name")
        end
    end
    
    # Get species information
    species_query = """
    SELECT 
        s.id, 
        s.name, 
        s.formula, 
        s.cas_number, 
        s.molecular_weight,
        s.metadata_json
    FROM 
        species s
    WHERE 
        s.id = ?
    """
    
    species_result = DuckDB.execute(conn, species_query, [species_id])
    species_df = DataFrame(species_result)
    
    # Get thermodynamic data sources for this species
    data_query = """
    SELECT 
        td.data_source,
        td.polynomial_type,
        td.temperature_min,
        td.temperature_max,
        td.reliability_score,
        ds.priority
    FROM 
        thermodynamic_data td
    JOIN 
        data_sources ds ON td.data_source = ds.name
    WHERE 
        td.species_id = ?
    ORDER BY 
        ds.priority DESC, 
        td.reliability_score DESC
    """
    
    data_result = DuckDB.execute(conn, data_query, [species_id])
    data_df = DataFrame(data_result)
    
    # Create result object
    result = Dict(
        "id" => species_df[1, :id],
        "name" => species_df[1, :name],
        "formula" => species_df[1, :formula],
        "cas_number" => species_df[1, :cas_number],
        "molecular_weight" => species_df[1, :molecular_weight]
    )
    
    # Parse metadata
    if !isempty(species_df[1, :metadata_json])
        try
            metadata = JSON.parse(species_df[1, :metadata_json])
            result["metadata"] = metadata
        catch e
            @warn "Failed to parse metadata JSON: $e"
        end
    end
    
    # Add data sources
    if size(data_df, 1) > 0
        sources = []
        
        for row in eachrow(data_df)
            source = Dict(
                "data_source" => row.data_source,
                "polynomial_type" => row.polynomial_type,
                "temperature_range" => [row.temperature_min, row.temperature_max],
                "reliability_score" => row.reliability_score,
                "priority" => row.priority
            )
            
            push!(sources, source)
        end
        
        result["data_sources"] = sources
    else
        result["data_sources"] = []
    end
    
    # Calculate coverage
    temp_range = config["general"]["temperature_range"]
    coverage = get_temperature_coverage(conn, species_name)
    
    if !isempty(coverage)
        # Calculate coverage percentage
        target_range = temp_range[2] - temp_range[1]
        covered_range = 0.0
        
        for row in eachrow(coverage)
            low = max(row.temperature_min, temp_range[1])
            high = min(row.temperature_max, temp_range[2])
            
            if high > low
                covered_range += high - low
            end
        end
        
        coverage_pct = min(covered_range / target_range * 100, 100.0)
        result["coverage_percentage"] = coverage_pct
        
        # Find gaps
        gaps = find_temperature_gaps(conn, species_name, temp_range)
        result["temperature_gaps"] = gaps
    else
        result["coverage_percentage"] = 0.0
        result["temperature_gaps"] = [temp_range]
    end
    
    return result
end

"""
    query_properties(conn::DuckDB.DB, species_name::String, temperature::Float64, config::Dict;
                    return_all_sources::Bool=false)

Query thermodynamic properties for a species at a specific temperature.
The data is progressively refined by traversing the hierarchy of data sources.
If return_all_sources is true, returns data from all sources instead of just the final refined result.
"""
function query_properties(conn::DuckDB.DB, species_name::String, temperature::Float64, config::Dict;
                         return_all_sources::Bool=false)
    # Check if in cache
    species_id = get_species_id(conn, species_name)
    
    if species_id === nothing
        error("Species not found: $species_name")
    end
    
    # If not requesting all sources, try to get from cache
    if !return_all_sources
        cached_data = get_from_cache(conn, species_id, temperature, "properties")
        
        if cached_data !== nothing
            return cached_data
        end
    end
    
    # Use the progressive refinement approach
    result, all_sources = progressively_refine_thermodynamic_data(conn, species_name, temperature, config)
    
    # Return all sources if requested
    if return_all_sources
        return Dict(
            "species_name" => species_name,
            "temperature" => temperature,
            "results" => all_sources
        )
    end
    
    return result
end

"""
    calculate_properties(conn::DuckDB.DB, species_name::String, temperature_range::Vector{Float64}, 
                       config::Dict; step::Float64=10.0)

Calculate thermodynamic properties over a temperature range.
"""
function calculate_properties(conn::DuckDB.DB, species_name::String, temperature_range::Vector{Float64}, 
                            config::Dict; step::Float64=10.0)
    # Validate temperature range
    if length(temperature_range) != 2 || temperature_range[1] >= temperature_range[2]
        error("Invalid temperature range: $temperature_range")
    end
    
    # Generate temperature points
    temperatures = temperature_range[1]:step:temperature_range[2]
    
    # Calculate properties at each temperature
    results = Dict[]
    
    for temp in temperatures
        try
            result = query_properties(conn, species_name, temp, config)
            push!(results, result)
        catch e
            @warn "Failed to calculate properties at $temp K: $e"
        end
    end
    
    # Extract property data for plotting
    temps = [result["temperature"] for result in results]
    
    cp_values = [result["properties"]["Cp"]["value"] for result in results]
    cp_uncertainties = [result["properties"]["Cp"]["uncertainty"] for result in results]
    
    h_values = [result["properties"]["H"]["value"] for result in results]
    h_uncertainties = [result["properties"]["H"]["uncertainty"] for result in results]
    
    s_values = [result["properties"]["S"]["value"] for result in results]
    s_uncertainties = [result["properties"]["S"]["uncertainty"] for result in results]
    
    g_values = [result["properties"]["G"]["value"] for result in results]
    g_uncertainties = [result["properties"]["G"]["uncertainty"] for result in results]
    
    # Create result structure
    return Dict(
        "species_name" => species_name,
        "temperature_range" => temperature_range,
        "temperatures" => temps,
        "properties" => Dict(
            "Cp" => Dict(
                "values" => cp_values,
                "uncertainties" => cp_uncertainties,
                "units" => "J/mol/K"
            ),
            "H" => Dict(
                "values" => h_values,
                "uncertainties" => h_uncertainties,
                "units" => "kJ/mol"
            ),
            "S" => Dict(
                "values" => s_values,
                "uncertainties" => s_uncertainties,
                "units" => "J/mol/K"
            ),
            "G" => Dict(
                "values" => g_values,
                "uncertainties" => g_uncertainties,
                "units" => "kJ/mol"
            )
        )
    )
end

"""
    list_available_species(conn::DuckDB.DB; limit::Int=100, search_term::String="")

List available species in the database.
"""
function list_available_species(conn::DuckDB.DB; limit::Int=100, search_term::String="")
    query = """
    SELECT 
        s.name,
        s.formula,
        s.cas_number,
        COUNT(td.id) AS data_count
    FROM 
        species s
    LEFT JOIN 
        thermodynamic_data td ON s.id = td.species_id
    """
    
    params = []
    
    if !isempty(search_term)
        query *= """
        WHERE 
            lower(s.name) LIKE lower(?) OR
            s.formula LIKE ? OR
            s.cas_number LIKE ?
        """
        
        search_pattern = "%$(search_term)%"
        push!(params, search_pattern, search_pattern, search_pattern)
    end
    
    query *= """
    GROUP BY 
        s.name, s.formula, s.cas_number
    ORDER BY 
        data_count DESC,
        s.name
    LIMIT ?
    """
    
    push!(params, limit)
    
    result = DuckDB.execute(conn, query, params)
    return DataFrame(result)
end

"""
    list_available_databases(conn::DuckDB.DB)

List all available thermodynamic databases.
"""
function list_available_databases(conn::DuckDB.DB)
    # Query data source information
    source_query = """
    SELECT 
        ds.name,
        ds.version,
        ds.url,
        ds.priority,
        ds.last_updated,
        COUNT(td.id) AS species_count
    FROM 
        data_sources ds
    LEFT JOIN 
        thermodynamic_data td ON ds.name = td.data_source
    GROUP BY 
        ds.name, ds.version, ds.url, ds.priority, ds.last_updated
    ORDER BY 
        ds.priority DESC
    """
    
    source_result = DuckDB.execute(conn, source_query)
    return DataFrame(source_result)
end