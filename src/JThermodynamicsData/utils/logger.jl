"""
Logging utilities for the JThermodynamicsData package.
"""

"""
    init_logger(level::String="info")

Initialize the logging system.
"""
function init_logger(level::String="info")
    # Parse log level
    log_level = parse_log_level(level)
    
    # Configure the default logger
    logger = ConsoleLogger(stderr, log_level)
    global_logger(logger)
    
    return logger
end

"""
    parse_log_level(level::String)

Parse a log level string into the corresponding Julia logging level.
"""
function parse_log_level(level::String)
    level_lower = lowercase(level)
    
    if level_lower == "debug"
        return Logging.Debug
    elseif level_lower == "info"
        return Logging.Info
    elseif level_lower == "warn"
        return Logging.Warn
    elseif level_lower == "error"
        return Logging.Error
    else
        @warn "Unknown log level: $level, using 'info'"
        return Logging.Info
    end
end

"""
    add_file_logger(file_path::String, level::String="info")

Add a file logger in addition to the console logger.
"""
function add_file_logger(file_path::String, level::String="info")
    # Ensure the directory exists
    mkpath(dirname(file_path))
    
    # Parse log level
    log_level = parse_log_level(level)
    
    # Create file logger
    file_logger = open(file_path, "a")
    
    # Create a TeeLogger that logs to both the console and the file
    tee_logger = TeeLogger(
        ConsoleLogger(stderr, log_level),
        SimpleLogger(file_logger, log_level)
    )
    
    global_logger(tee_logger)
    
    return tee_logger, file_logger
end

"""
    close_file_logger(file_logger::IO)

Close a file logger opened with add_file_logger.
"""
function close_file_logger(file_logger::IO)
    close(file_logger)
    
    # Reset to console logger
    init_logger()
    
    return nothing
end

"""
    log_operation(operation::String, params::Dict)

Log an operation with its parameters.
"""
function log_operation(operation::String, params::Dict)
    # Format parameters
    params_str = join(["$k=$(repr(v))" for (k, v) in params], ", ")
    
    @info "Operation: $operation($params_str)"
    
    return nothing
end

"""
    log_database_query(conn::DuckDB.DB, query::String, params::Vector)

Log a database query.
"""
function log_database_query(conn::DuckDB.DB, query::String, params::Vector)
    # Only log in debug mode
    @debug "SQL Query: $query" params=params
    
    return nothing
end

"""
    log_error(error_msg::String, context::Dict=Dict())

Log an error with context information.
"""
function log_error(error_msg::String, context::Dict=Dict())
    # Format context
    context_str = isempty(context) ? "" : " ($(join(["$k=$(repr(v))" for (k, v) in context], ", ")))"
    
    @error "Error: $error_msg$context_str"
    
    return nothing
end