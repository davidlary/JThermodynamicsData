"""
Logging utilities for the JThermodynamicsData package.

This module provides comprehensive logging capabilities for tracking operations,
timing performance, recording pipeline stages, and capturing errors throughout
the thermodynamic data processing pipeline.
"""

using Dates
using Printf

# Store global timing benchmarks
const _timing_benchmarks = Dict{String, Dict{String, Any}}()

"""
    init_logger(level::String="info", log_file::Union{Nothing,String}=nothing)

Initialize the logging system with optional file output.

# Arguments
- `level::String`: The log level to use ("debug", "info", "warn", "error")
- `log_file::Union{Nothing,String}`: Optional file path for logging output

# Returns
- The configured logger
"""
function init_logger(level::String="info", log_file::Union{Nothing,String}=nothing)
    # Parse log level
    log_level = parse_log_level(level)
    
    # Configure the default logger based on whether we have a log file
    if log_file === nothing
        logger = ConsoleLogger(stderr, log_level)
        global_logger(logger)
        return logger
    else
        return add_file_logger(log_file, level)
    end
end

"""
    parse_log_level(level::String)

Parse a log level string into the corresponding Julia logging level.

# Arguments
- `level::String`: The log level string ("debug", "info", "warn", "error")

# Returns
- The corresponding Julia Logging level
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

# Arguments
- `file_path::String`: The path to write logs to
- `level::String`: The log level to use

# Returns
- Tuple of (tee_logger, file_logger) where file_logger should be closed when done
"""
function add_file_logger(file_path::String, level::String="info")
    # Ensure the directory exists
    mkpath(dirname(file_path))
    
    # Parse log level
    log_level = parse_log_level(level)
    
    # Create file logger with write mode - this will truncate if exists
    file_logger = open(file_path, "w")
    
    # Write header to log file
    println(file_logger, "="^80)
    println(file_logger, "JThermodynamicsData Log - $(now())")
    println(file_logger, "Log Level: $(uppercase(level))")
    println(file_logger, "="^80)
    println(file_logger)
    
    # Ensure header is written immediately
    flush(file_logger)
    
    # Create a custom LoggingExtras formatter to ensure flushing
    format_logger = TransformerLogger(SimpleLogger(file_logger, log_level)) do log
        # Flush the file after each log message
        flush(file_logger)
        return log
    end
    
    # Create a TeeLogger that logs to both the console and the file
    tee_logger = TeeLogger(
        ConsoleLogger(stderr, log_level),
        format_logger
    )
    
    global_logger(tee_logger)
    
    return tee_logger, file_logger
end

"""
    close_file_logger(file_logger::IO)

Close a file logger opened with add_file_logger.

# Arguments
- `file_logger::IO`: The file logger to close
"""
function close_file_logger(file_logger::IO)
    # Write footer to log file
    println(file_logger)
    println(file_logger, "="^80)
    println(file_logger, "Log Closed - $(now())")
    println(file_logger, "="^80)
    
    # Ensure footer is written
    flush(file_logger)
    
    # Close file
    close(file_logger)
    
    # Reset to console logger
    init_logger()
    
    return nothing
end

"""
    log_pipeline_start(pipeline_name::String, config::Dict)

Log the start of a pipeline execution with configuration details.

# Arguments
- `pipeline_name::String`: The name of the pipeline
- `config::Dict`: Configuration parameters for the pipeline
"""
function log_pipeline_start(pipeline_name::String, config::Dict)
    @info "======== PIPELINE START: $pipeline_name ========"
    
    # Log key configuration parameters if available
    # Check if this is a full project config or a simple config dict
    if haskey(config, "general") && isa(config["general"], Dict)
        # This is a full project config
        @info "Pipeline Configuration:" pipeline_name config=Dict(
            "temperature_range" => get(config["general"], "temperature_range", "N/A"),
            "database_path" => get(config["general"], "database_path", "N/A"),
            "log_level" => get(config["general"], "log_level", "info"),
            "enabled_sources" => haskey(config, "data_sources") ? 
                [src["name"] for src in config["data_sources"] if get(src, "enabled", false)] : 
                []
        )
    else
        # This is a simple config dict, just log it as-is
        @info "Pipeline Configuration:" pipeline_name config=config
    end
    
    # Start a global timer for the pipeline
    global _timing_benchmarks
    _timing_benchmarks[pipeline_name] = Dict(
        "start_time" => now(),
        "stages" => Dict{String, Dict{String, Any}}()
    )
    
    return nothing
end

"""
    log_pipeline_end(pipeline_name::String, results::Dict=Dict())

Log the end of a pipeline execution with timing information.

# Arguments
- `pipeline_name::String`: The name of the pipeline
- `results::Dict`: Optional results summary
"""
function log_pipeline_end(pipeline_name::String, results::Dict=Dict())
    if !haskey(_timing_benchmarks, pipeline_name)
        @warn "No start time found for pipeline: $pipeline_name"
        return nothing
    end
    
    # Calculate elapsed time
    global _timing_benchmarks
    start_time = _timing_benchmarks[pipeline_name]["start_time"]
    elapsed = now() - start_time
    
    # Log completion message with timing
    @info "======== PIPELINE COMPLETE: $pipeline_name ========"
    @info "Pipeline Execution Time: $(format_time_duration(elapsed))"
    
    # Log summary of results if provided
    if !isempty(results)
        @info "Pipeline Results Summary:" results
    end
    
    # Log timing for each stage
    if haskey(_timing_benchmarks[pipeline_name], "stages") && !isempty(_timing_benchmarks[pipeline_name]["stages"])
        stages = _timing_benchmarks[pipeline_name]["stages"]
        @info "Stage Timing Summary:"
        for (stage_name, stage_data) in sort(collect(stages), by=x->x[2]["start_time"])
            if haskey(stage_data, "end_time")
                stage_elapsed = stage_data["end_time"] - stage_data["start_time"]
                @info "  - $stage_name: $(format_time_duration(stage_elapsed))"
            end
        end
    end
    
    return nothing
end

"""
    log_stage_start(pipeline_name::String, stage_name::String, params::Dict=Dict())

Log the start of a pipeline stage with parameters.

# Arguments
- `pipeline_name::String`: The name of the pipeline
- `stage_name::String`: The name of the stage
- `params::Dict`: Parameters for the stage
"""
function log_stage_start(pipeline_name::String, stage_name::String, params::Dict=Dict())
    global _timing_benchmarks
    
    # Create pipeline entry if it doesn't exist
    if !haskey(_timing_benchmarks, pipeline_name)
        _timing_benchmarks[pipeline_name] = Dict(
            "start_time" => now(),
            "stages" => Dict{String, Dict{String, Any}}()
        )
    end
    
    # Record stage start time
    if !haskey(_timing_benchmarks[pipeline_name], "stages")
        _timing_benchmarks[pipeline_name]["stages"] = Dict{String, Dict{String, Any}}()
    end
    
    _timing_benchmarks[pipeline_name]["stages"][stage_name] = Dict(
        "start_time" => now(),
        "params" => params
    )
    
    # Format parameters
    params_str = isempty(params) ? "" : join([" $k=$(repr(v))" for (k, v) in params], ",")
    
    @info "Stage Start: $stage_name$params_str"
    
    return nothing
end

"""
    log_stage_end(pipeline_name::String, stage_name::String, results::Dict=Dict())

Log the end of a pipeline stage with results.

# Arguments
- `pipeline_name::String`: The name of the pipeline
- `stage_name::String`: The name of the stage
- `results::Dict`: Results of the stage execution
"""
function log_stage_end(pipeline_name::String, stage_name::String, results::Dict=Dict())
    global _timing_benchmarks
    
    # Verify we have timing data for this stage
    if !haskey(_timing_benchmarks, pipeline_name) || 
       !haskey(_timing_benchmarks[pipeline_name], "stages") ||
       !haskey(_timing_benchmarks[pipeline_name]["stages"], stage_name)
        @warn "No start time found for stage: $stage_name in pipeline: $pipeline_name"
        return nothing
    end
    
    # Record stage end time
    stage_data = _timing_benchmarks[pipeline_name]["stages"][stage_name]
    stage_data["end_time"] = now()
    
    # Calculate elapsed time
    elapsed = stage_data["end_time"] - stage_data["start_time"]
    
    # Format results
    results_str = isempty(results) ? "" : join([" $k=$(repr(v))" for (k, v) in results], ",")
    
    @info "Stage Complete: $stage_name ($(format_time_duration(elapsed)))$results_str"
    
    return nothing
end

"""
    log_operation(operation::String, params::Dict)

Log an operation with its parameters.

# Arguments
- `operation::String`: The name of the operation
- `params::Dict`: Parameters for the operation
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

# Arguments
- `conn::DuckDB.DB`: The database connection
- `query::String`: The SQL query
- `params::Vector`: Query parameters
"""
function log_database_query(conn::DuckDB.DB, query::String, params::Vector)
    # Only log in debug mode
    @debug "SQL Query: $query" params=params
    
    return nothing
end

"""
    log_error(error_msg::String, context::Dict=Dict())

Log an error with context information.

# Arguments
- `error_msg::String`: The error message
- `context::Dict`: Context information about the error
"""
function log_error(error_msg::String, context::Dict=Dict())
    # Format context
    context_str = isempty(context) ? "" : " ($(join(["$k=$(repr(v))" for (k, v) in context], ", ")))"
    
    @error "Error: $error_msg$context_str"
    
    return nothing
end

"""
    log_species_processing(species_name::String, temperature_range::Vector{<:Real}, step::Real)

Log the start of processing for a specific species.

# Arguments
- `species_name::String`: The name of the species
- `temperature_range::Vector{<:Real}`: The temperature range being used
- `step::Real`: The temperature step size
"""
function log_species_processing(species_name::String, temperature_range::Vector{<:Real}, step::Real)
    @info "Processing Species: $species_name" temperature_range=temperature_range step=step
    return nothing
end

"""
    log_parser_operation(parser_name::String, source_name::String, file_path::String)

Log a parser operation for a specific data source.

# Arguments
- `parser_name::String`: The name of the parser
- `source_name::String`: The name of the data source
- `file_path::String`: The path to the file being parsed
"""
function log_parser_operation(parser_name::String, source_name::String, file_path::String)
    @info "Running Parser: $parser_name" source=source_name file=file_path
    return nothing
end

"""
    log_timing_benchmark(category::String, operation::String, elapsed::TimePeriod)

Log a timing benchmark for performance analysis.

# Arguments
- `category::String`: The category of the operation
- `operation::String`: The specific operation
- `elapsed::TimePeriod`: The elapsed time
"""
function log_timing_benchmark(category::String, operation::String, elapsed::TimePeriod)
    @info "Performance: [$category] $operation took $(format_time_duration(elapsed))"
    return nothing
end

"""
    format_time_duration(elapsed::TimePeriod)

Format a time duration for display in logs.

# Arguments
- `elapsed::TimePeriod`: The elapsed time period

# Returns
- A formatted string representation of the duration
"""
function format_time_duration(elapsed::TimePeriod)
    total_milliseconds = Dates.value(convert(Millisecond, elapsed))
    
    if total_milliseconds < 1000
        return "$(total_milliseconds)ms"
    elseif total_milliseconds < 60000
        return @sprintf("%.2fs", total_milliseconds / 1000)
    elseif total_milliseconds < 3600000
        minutes = div(total_milliseconds, 60000)
        seconds = (total_milliseconds % 60000) / 1000
        return @sprintf("%dm %.2fs", minutes, seconds)
    else
        hours = div(total_milliseconds, 3600000)
        minutes = div(total_milliseconds % 3600000, 60000)
        seconds = (total_milliseconds % 60000) / 1000
        return @sprintf("%dh %dm %.2fs", hours, minutes, seconds)
    end
end

"""
    with_timing(f::Function, category::String, operation::String)

Execute a function and log its execution time.

# Arguments
- `f::Function`: The function to execute
- `category::String`: The category of the operation
- `operation::String`: The specific operation

# Returns
- The result of the function call
"""
function with_timing(f::Function, category::String, operation::String)
    start_time = now()
    result = f()
    elapsed = now() - start_time
    
    log_timing_benchmark(category, operation, elapsed)
    
    return result
end

"""
    create_log_context()

Create a new logging context with a unique ID for tracking related log entries.

# Returns
- A Dict with context information including a UUID
"""
function create_log_context()
    context_id = rand(UInt64)
    return Dict(
        "context_id" => "ctx-$(string(context_id, base=16))",
        "created_at" => now()
    )
end