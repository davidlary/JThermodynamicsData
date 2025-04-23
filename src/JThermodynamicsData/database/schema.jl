"""
Database schema definitions and initialization functions.
"""

"""
    create_schema(conn::DuckDB.DB)

Create the database schema for thermodynamic data.
"""
function create_schema(conn::DuckDB.DB)
    # Create sequences for auto-increment IDs
    DuckDB.execute(conn, "CREATE SEQUENCE IF NOT EXISTS species_id_seq")
    DuckDB.execute(conn, "CREATE SEQUENCE IF NOT EXISTS thermodynamic_data_id_seq")
    DuckDB.execute(conn, "CREATE SEQUENCE IF NOT EXISTS data_sources_id_seq")
    DuckDB.execute(conn, "CREATE SEQUENCE IF NOT EXISTS data_cache_id_seq")
    DuckDB.execute(conn, "CREATE SEQUENCE IF NOT EXISTS database_version_id_seq")
    
    # Species table
    DuckDB.execute(conn, """
    CREATE TABLE IF NOT EXISTS species (
        id INTEGER PRIMARY KEY,
        name VARCHAR NOT NULL,
        formula VARCHAR NOT NULL,
        cas_number VARCHAR,
        molecular_weight DOUBLE,
        inchi VARCHAR,
        smiles VARCHAR,
        metadata_json VARCHAR,
        CONSTRAINT species_name_unique UNIQUE (name)
    )
    """)
    
    # Create index on species name (case insensitive)
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_species_name ON species (lower(name))")
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_species_formula ON species (formula)")
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_species_cas ON species (cas_number)")
    
    # Thermodynamic data table
    DuckDB.execute(conn, """
    CREATE TABLE IF NOT EXISTS thermodynamic_data (
        id INTEGER PRIMARY KEY,
        species_id INTEGER NOT NULL,
        data_source VARCHAR NOT NULL,
        polynomial_type VARCHAR NOT NULL,
        temperature_min DOUBLE NOT NULL,
        temperature_max DOUBLE NOT NULL,
        reliability_score DOUBLE DEFAULT 0.0,
        data_json VARCHAR NOT NULL,
        uncertainty_json VARCHAR,
        date_added TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        date_modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (species_id) REFERENCES species(id),
        CONSTRAINT thermo_data_unique UNIQUE (species_id, data_source, polynomial_type)
    )
    """)
    
    # Create indexes on thermodynamic_data table
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_thermo_species_id ON thermodynamic_data (species_id)")
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_thermo_source ON thermodynamic_data (data_source)")
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_thermo_reliability ON thermodynamic_data (reliability_score DESC)")
    DuckDB.execute(conn, "CREATE INDEX IF NOT EXISTS idx_thermo_temp_range ON thermodynamic_data (temperature_min, temperature_max)")
    
    # Data source metadata table
    DuckDB.execute(conn, """
    CREATE TABLE IF NOT EXISTS data_sources (
        id INTEGER PRIMARY KEY,
        name VARCHAR NOT NULL,
        version VARCHAR,
        url VARCHAR,
        description VARCHAR,
        priority INTEGER NOT NULL,
        last_updated TIMESTAMP,
        metadata_json VARCHAR,
        CONSTRAINT source_name_unique UNIQUE (name)
    )
    """)
    
    # Version tracking table
    DuckDB.execute(conn, """
    CREATE TABLE IF NOT EXISTS database_version (
        id INTEGER PRIMARY KEY,
        version INTEGER NOT NULL,
        schema_version INTEGER NOT NULL,
        description VARCHAR,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """)
    
    # Check if version record exists, then insert if it doesn't
    result = DuckDB.execute(conn, "SELECT COUNT(*) as count FROM database_version")
    df = DataFrame(result)
    
    if size(df, 1) == 0 || df[1, :count] == 0
        # Insert initial version record
        DuckDB.execute(conn, """
        INSERT INTO database_version (id, version, schema_version, description)
        VALUES (nextval('database_version_id_seq'), 1, 1, 'Initial database creation')
        """)
    end
    
    # Cache table for frequently accessed data
    DuckDB.execute(conn, """
    CREATE TABLE IF NOT EXISTS data_cache (
        id INTEGER PRIMARY KEY,
        species_id INTEGER NOT NULL,
        temperature_min DOUBLE NOT NULL,
        temperature_max DOUBLE NOT NULL,
        cache_type VARCHAR NOT NULL,
        cache_data VARCHAR NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (species_id) REFERENCES species(id),
        CONSTRAINT cache_unique UNIQUE (species_id, temperature_min, temperature_max, cache_type)
    )
    """)
    
    # In DuckDB, we'll skip the function for now
    # Function definition varies by version, and we can calculate this in Julia instead
    
    # View for combined species and thermodynamic data
    DuckDB.execute(conn, """
    CREATE OR REPLACE VIEW v_thermodynamic_data AS
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
        td.date_modified
    FROM
        species s
    JOIN
        thermodynamic_data td ON s.id = td.species_id
    """)
    
    # Initialize data sources from config
    init_data_sources(conn)
    
    return true
end

"""
    init_data_sources(conn::DuckDB.DB)

Initialize data sources table from configuration.
"""
function init_data_sources(conn::DuckDB.DB; config=nothing)
    if config === nothing
        config_path = joinpath(dirname(dirname(dirname(dirname(@__FILE__)))), "config", "settings.yaml")
        config = load_config(config_path)
    end
    
    # Clear existing sources
    DuckDB.execute(conn, "DELETE FROM data_sources")
    
    # Insert data sources from config
    for (i, source) in enumerate(config["data_sources"])
        name = source["name"]
        priority = source["priority"]
        url = get(source, "url", "")
        
        # Use direct execution instead of prepare statement
        metadata = Dict(
            "format" => source["format"],
            "enabled" => source["enabled"],
            "reliability_score" => get(source, "reliability_score", 0.0),
            "refresh_interval_days" => get(source, "refresh_interval_days", 90)
        )
        
        metadata_json = JSON.json(metadata)
        
        DuckDB.execute(conn, """
        INSERT INTO data_sources (id, name, version, url, priority, metadata_json)
        VALUES (?, ?, '1.0', ?, ?, ?)
        """, [i, name, url, priority, metadata_json])
    end
    
    return true
end

"""
    get_schema_version(conn::DuckDB.DB)

Get the current schema version.
"""
function get_schema_version(conn::DuckDB.DB)
    result = DuckDB.execute(conn, """
    SELECT schema_version
    FROM database_version
    ORDER BY id DESC
    LIMIT 1
    """)
    
    df = DataFrame(result)
    if size(df, 1) > 0
        return df[1, :schema_version]
    else
        return 0
    end
end

"""
    upgrade_schema(conn::DuckDB.DB, target_version::Int)

Upgrade the schema to the target version.
"""
function upgrade_schema(conn::DuckDB.DB, target_version::Int)
    current_version = get_schema_version(conn)
    
    if current_version >= target_version
        return true
    end
    
    # For future schema upgrades, implement version-specific migrations here
    if current_version < 2 && target_version >= 2
        # Example: Upgrade from version 1 to 2
        # DuckDB.execute(conn, "ALTER TABLE species ADD COLUMN new_column VARCHAR")
    end
    
    # Update version record
    DuckDB.execute(conn, """
    INSERT INTO database_version (id, version, schema_version, description)
    VALUES (nextval('database_version_id_seq'), ?, ?, ?)
    """, [current_version + 1, target_version, "Upgraded schema to version $target_version"])
    
    return true
end