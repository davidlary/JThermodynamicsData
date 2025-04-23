"""
Database connection management functions.
"""

"""
    connect_database(db_path::String)

Connect to an existing database without reinitializing it.
"""
function connect_database(db_path::String)
    if !isfile(db_path)
        error("Database file not found: $db_path")
    end
    
    # Connect to existing DuckDB database
    conn = DuckDB.DB(db_path)
    
    return conn
end

"""
    init_database(config::Dict)

Initialize the database connection and schema.
"""
function init_database(config::Dict)
    db_path = config["general"]["database_path"]
    
    # Create parent directory if it doesn't exist
    mkpath(dirname(db_path))
    
    # Remove existing database if it exists to avoid constraint errors
    if isfile(db_path)
        rm(db_path)
    end
    
    # Connect to DuckDB database
    conn = DuckDB.DB(db_path)
    
    # Create schema if needed
    create_schema(conn)
    
    # Set database configuration
    set_database_config(conn, config["database"])
    
    return conn
end

"""
    set_database_config(conn::DuckDB.DB, db_config::Dict)

Set database configuration parameters.
"""
function set_database_config(conn::DuckDB.DB, db_config::Dict)
    # Set cache size
    cache_size = get(db_config, "cache_size_mb", 1024)
    DuckDB.execute(conn, "PRAGMA memory_limit='$(cache_size)MB'")
    
    # Set query timeout
    timeout = get(db_config, "query_timeout_seconds", 30)
    # Skip parallel query setting - varies by DuckDB version
    
    return true
end

"""
    close_database(conn::DuckDB.DB, vacuum::Bool=true)

Close the database connection, optionally vacuuming the database.
"""
function close_database(conn::DuckDB.DB, vacuum::Bool=true)
    if vacuum
        DuckDB.execute(conn, "VACUUM")
    end
    
    # DuckDB.disconnect is no longer needed in newer versions
    # The connection will be cleaned up by the garbage collector
    return true
end

"""
    transaction(func::Function, conn::DuckDB.DB)

Execute a function within a database transaction.
"""
function transaction(func::Function, conn::DuckDB.DB)
    DuckDB.execute(conn, "BEGIN TRANSACTION")
    
    try
        result = func()
        DuckDB.execute(conn, "COMMIT")
        return result
    catch e
        DuckDB.execute(conn, "ROLLBACK")
        rethrow(e)
    end
end

"""
    vacuum_database(conn::DuckDB.DB)

Vacuum the database to optimize storage.
"""
function vacuum_database(conn::DuckDB.DB)
    DuckDB.execute(conn, "VACUUM")
    return true
end

"""
    analyze_database(conn::DuckDB.DB)

Analyze the database to update statistics.
"""
function analyze_database(conn::DuckDB.DB)
    DuckDB.execute(conn, "ANALYZE")
    return true
end

"""
    backup_database(conn::DuckDB.DB, backup_path::String)

Backup the database to the specified path.
"""
function backup_database(conn::DuckDB.DB, backup_path::String)
    # Create parent directory if it doesn't exist
    mkpath(dirname(backup_path))
    
    # Create a backup
    DuckDB.execute(conn, "EXPORT DATABASE '$(backup_path)'")
    
    return true
end

"""
    restore_database(conn::DuckDB.DB, backup_path::String)

Restore the database from the specified backup.
"""
function restore_database(conn::DuckDB.DB, backup_path::String)
    if !isfile(backup_path)
        error("Backup file not found: $backup_path")
    end
    
    # Restore from backup
    DuckDB.execute(conn, "IMPORT DATABASE '$(backup_path)'")
    
    return true
end