"""
Database query functions for thermodynamic data.
"""

"""
    get_species_id(conn::DuckDB.DB, species_name::String)

Get the ID for a species by name.
"""
function get_species_id(conn::DuckDB.DB, species_name::String)
    result = DuckDB.execute(conn, """
    SELECT id FROM species
    WHERE lower(name) = lower(?)
    LIMIT 1
    """, [species_name])
    
    df = DataFrame(result)
    if size(df, 1) > 0
        return df[1, :id]
    else
        return nothing
    end
end

"""
    get_species_by_formula(conn::DuckDB.DB, formula::String)

Get all species with the given formula.
"""
function get_species_by_formula(conn::DuckDB.DB, formula::String)
    result = DuckDB.execute(conn, """
    SELECT id, name, formula, cas_number, molecular_weight
    FROM species
    WHERE formula = ?
    """, [formula])
    
    return DataFrame(result)
end

"""
    get_species_by_cas(conn::DuckDB.DB, cas_number::String)

Get a species by CAS registry number.
"""
function get_species_by_cas(conn::DuckDB.DB, cas_number::String)
    result = DuckDB.execute(conn, """
    SELECT id, name, formula, cas_number, molecular_weight
    FROM species
    WHERE cas_number = ?
    LIMIT 1
    """, [cas_number])
    
    df = DataFrame(result)
    if size(df, 1) > 0
        return df[1, :]
    else
        return nothing
    end
end

"""
    add_species(conn::DuckDB.DB, name::String, formula::String, cas_number::String="", 
                molecular_weight::Float64=0.0, metadata::Dict=Dict())

Add a new species to the database.
"""
function add_species(conn::DuckDB.DB, name::String, formula::String, cas_number::String="", 
                     molecular_weight::Float64=0.0, metadata::Dict=Dict())
    # First check if species already exists
    species_id = get_species_id(conn, name)
    if species_id !== nothing
        return species_id
    end
    
    # Get max ID to generate a new ID
    max_id_result = DuckDB.execute(conn, "SELECT COALESCE(MAX(id), 0) as max_id FROM species")
    max_id_df = DataFrame(max_id_result)
    new_id = max_id_df[1, :max_id] + 1
    
    # Insert new species with explicit ID
    metadata_json = isempty(metadata) ? "{}" : JSON.json(metadata)
    
    # Use parameterized query directly without prepare
    query = """
    INSERT INTO species (id, name, formula, cas_number, molecular_weight, metadata_json)
    VALUES (?, ?, ?, ?, ?, ?)
    RETURNING id
    """
    
    result = DuckDB.execute(conn, query, [new_id, name, formula, cas_number, molecular_weight, metadata_json])
    df = DataFrame(result)
    
    return df[1, :id]
end

"""
    add_thermodynamic_data(conn::DuckDB.DB, species_id::Int, data_source::String,
                          polynomial_type::String, temp_min::Float64, temp_max::Float64,
                          data::Dict, uncertainty::Dict=Dict(), reliability_score::Float64=0.0)

Add thermodynamic data for a species.
"""
function add_thermodynamic_data(conn::DuckDB.DB, species_id::Int, data_source::String,
                               polynomial_type::String, temp_min::Float64, temp_max::Float64,
                               data::Dict, uncertainty::Dict=Dict(), reliability_score::Float64=0.0)
    # Convert data to JSON
    data_json = JSON.json(data)
    uncertainty_json = isempty(uncertainty) ? "{}" : JSON.json(uncertainty)
    
    # Check if data already exists
    result = DuckDB.execute(conn, """
    SELECT id FROM thermodynamic_data
    WHERE species_id = ? AND data_source = ? AND polynomial_type = ?
    LIMIT 1
    """, [species_id, data_source, polynomial_type])
    
    df = DataFrame(result)
    
    if size(df, 1) > 0
        # Update existing record
        record_id = df[1, :id]
        
        # Use parameterized query directly
        query = """
        UPDATE thermodynamic_data
        SET temperature_min = ?, temperature_max = ?, reliability_score = ?,
            data_json = ?, uncertainty_json = ?, date_modified = CURRENT_TIMESTAMP
        WHERE id = ?
        """
        
        DuckDB.execute(conn, query, [temp_min, temp_max, reliability_score, data_json, uncertainty_json, record_id])
        return record_id
    else
        # Get max ID to generate a new ID
        max_id_result = DuckDB.execute(conn, "SELECT COALESCE(MAX(id), 0) as max_id FROM thermodynamic_data")
        max_id_df = DataFrame(max_id_result)
        new_id = max_id_df[1, :max_id] + 1
        
        # Insert new record with explicit ID
        query = """
        INSERT INTO thermodynamic_data 
        (id, species_id, data_source, polynomial_type, temperature_min, temperature_max, 
         reliability_score, data_json, uncertainty_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        RETURNING id
        """
        
        result = DuckDB.execute(conn, query, [
            new_id, species_id, data_source, polynomial_type, temp_min, temp_max,
            reliability_score, data_json, uncertainty_json
        ])
        
        df = DataFrame(result)
        return df[1, :id]
    end
end

"""
    get_thermodynamic_data(conn::DuckDB.DB, species_name::String; 
                         temp_min::Float64=0.0, temp_max::Float64=1e6,
                         data_source::Union{String, Nothing}=nothing)

Get thermodynamic data for a species within a temperature range.
"""
function get_thermodynamic_data(conn::DuckDB.DB, species_name::String; 
                              temp_min::Float64=0.0, temp_max::Float64=1e6,
                              data_source::Union{String, Nothing}=nothing)
    # Build query conditions
    conditions = ["s.name = ?"]
    params = [species_name]
    
    if temp_min > 0
        push!(conditions, "td.temperature_max >= ?")
        push!(params, Float64(temp_min))
    end
    
    if temp_max < 1e6
        push!(conditions, "td.temperature_min <= ?")
        push!(params, Float64(temp_max))
    end
    
    if data_source !== nothing
        push!(conditions, "td.data_source = ?")
        push!(params, data_source)
    end
    
    # Join conditions with AND
    where_clause = join(conditions, " AND ")
    
    # Execute query
    query = """
    SELECT 
        s.name AS species_name,
        s.formula,
        s.cas_number,
        s.molecular_weight,
        td.data_source,
        td.polynomial_type,
        td.temperature_min,
        td.temperature_max,
        td.reliability_score,
        td.data_json,
        td.uncertainty_json,
        ds.priority,
        td.date_modified
    FROM 
        species s
    JOIN 
        thermodynamic_data td ON s.id = td.species_id
    JOIN 
        data_sources ds ON td.data_source = ds.name
    WHERE 
        $where_clause
    ORDER BY 
        ds.priority DESC, 
        td.reliability_score DESC,
        td.date_modified DESC
    """
    
    result = DuckDB.execute(conn, query, params)
    return DataFrame(result)
end

"""
    get_best_thermodynamic_data(conn::DuckDB.DB, species_name::String, temperature::Float64)

Get the best thermodynamic data for a species at a specific temperature.
"""
function get_best_thermodynamic_data(conn::DuckDB.DB, species_name::String, temperature::Float64)
    df = get_thermodynamic_data(conn, species_name; temp_min=temperature, temp_max=temperature)
    
    if size(df, 1) > 0
        # Get first row (highest priority/reliability)
        row = df[1, :]
        
        # Parse JSON data
        data = JSON.parse(row.data_json)
        uncertainty = JSON.parse(row.uncertainty_json)
        
        return Dict(
            "species_name" => row.species_name,
            "formula" => row.formula,
            "cas_number" => row.cas_number,
            "molecular_weight" => row.molecular_weight,
            "data_source" => row.data_source,
            "polynomial_type" => row.polynomial_type,
            "temperature_range" => [row.temperature_min, row.temperature_max],
            "reliability_score" => row.reliability_score,
            "data" => data,
            "uncertainty" => uncertainty,
            "date_modified" => row.date_modified
        )
    else
        return nothing
    end
end

"""
    search_species(conn::DuckDB.DB, search_term::String; limit::Int=100)

Search for species by name, formula, or CAS number.
"""
function search_species(conn::DuckDB.DB, search_term::String; limit::Int=100)
    search_pattern = "%$(search_term)%"
    
    result = DuckDB.execute(conn, """
    SELECT 
        s.id, 
        s.name, 
        s.formula, 
        s.cas_number, 
        s.molecular_weight,
        COUNT(td.id) AS data_count
    FROM 
        species s
    LEFT JOIN 
        thermodynamic_data td ON s.id = td.species_id
    WHERE 
        lower(s.name) LIKE lower(?) OR
        s.formula LIKE ? OR
        s.cas_number LIKE ?
    GROUP BY 
        s.id, s.name, s.formula, s.cas_number, s.molecular_weight
    ORDER BY 
        data_count DESC,
        CASE WHEN lower(s.name) = lower(?) THEN 1
             WHEN lower(s.name) LIKE lower(?) THEN 2
             ELSE 3
        END
    LIMIT ?
    """, [search_pattern, search_pattern, search_pattern, search_term, search_pattern, limit])
    
    return DataFrame(result)
end

"""
    list_data_sources(conn::DuckDB.DB)

List all available data sources with their metadata.
"""
function list_data_sources(conn::DuckDB.DB)
    result = DuckDB.execute(conn, """
    SELECT 
        name, 
        version, 
        url, 
        priority,
        last_updated,
        metadata_json
    FROM 
        data_sources
    ORDER BY 
        priority DESC
    """)
    
    df = DataFrame(result)
    
    # Parse metadata JSON
    if size(df, 1) > 0
        df.metadata = map(row -> JSON.parse(row), df.metadata_json)
        select!(df, Not(:metadata_json))
    end
    
    return df
end

"""
    update_data_source_metadata(conn::DuckDB.DB, source_name::String, metadata::Dict)

Update metadata for a data source.
"""
function update_data_source_metadata(conn::DuckDB.DB, source_name::String, metadata::Dict)
    metadata_json = JSON.json(metadata)
    
    query = """
    UPDATE data_sources
    SET metadata_json = ?, last_updated = CURRENT_TIMESTAMP
    WHERE name = ?
    """
    
    DuckDB.execute(conn, query, [metadata_json, source_name])
    
    return true
end

"""
    update_data_source_version(conn::DuckDB.DB, source_name::String, version::String)

Update version for a data source.
"""
function update_data_source_version(conn::DuckDB.DB, source_name::String, version::String)
    query = """
    UPDATE data_sources
    SET version = ?, last_updated = CURRENT_TIMESTAMP
    WHERE name = ?
    """
    
    DuckDB.execute(conn, query, [version, source_name])
    
    return true
end

"""
    delete_thermodynamic_data(conn::DuckDB.DB, species_name::String, data_source::String)

Delete thermodynamic data for a species from a specific source.
"""
function delete_thermodynamic_data(conn::DuckDB.DB, species_name::String, data_source::String)
    species_id = get_species_id(conn, species_name)
    
    if species_id === nothing
        return false
    end
    
    query = """
    DELETE FROM thermodynamic_data
    WHERE species_id = ? AND data_source = ?
    """
    
    DuckDB.execute(conn, query, [species_id, data_source])
    
    return true
end

"""
    get_temperature_coverage(conn::DuckDB.DB, species_name::String)

Get temperature coverage for a species across all data sources.
"""
function get_temperature_coverage(conn::DuckDB.DB, species_name::String)
    species_id = get_species_id(conn, species_name)
    
    if species_id === nothing
        return []
    end
    
    result = DuckDB.execute(conn, """
    SELECT 
        data_source,
        temperature_min,
        temperature_max
    FROM 
        thermodynamic_data
    WHERE 
        species_id = ?
    ORDER BY 
        temperature_min
    """, [species_id])
    
    return DataFrame(result)
end

"""
    find_temperature_gaps(conn::DuckDB.DB, species_name::String, target_range::Vector{Float64})

Find gaps in temperature coverage for a species.
"""
function find_temperature_gaps(conn::DuckDB.DB, species_name::String, target_range::Vector{Float64})
    coverage = get_temperature_coverage(conn, species_name)
    
    if isempty(coverage)
        return [target_range]
    end
    
    # Sort by temperature_min
    sort!(coverage, :temperature_min)
    
    gaps = []
    current_max = target_range[1]
    
    for i in 1:size(coverage, 1)
        row = coverage[i, :]
        
        if row.temperature_min > current_max
            # Found a gap
            push!(gaps, [current_max, row.temperature_min])
        end
        
        current_max = max(current_max, row.temperature_max)
    end
    
    if current_max < target_range[2]
        # Gap at the end
        push!(gaps, [current_max, target_range[2]])
    end
    
    return gaps
end