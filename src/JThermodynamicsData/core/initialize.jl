"""
    Core initialization module for JThermodynamicsData

This module provides functions to initialize the entire system,
including configuration loading, database connection, and setting up
the hierarchical source traversal for thermodynamic data.
"""

"""
    initialize(config_path::String="config/settings.yaml")

Initialize the JThermodynamicsData package with configuration.
Loads all configuration settings and connects to the database.

# Arguments
- `config_path::String`: Path to the main configuration file

# Returns
- A tuple of (config, database_connection)
"""
function initialize(config_path::String="config/settings.yaml")
    # Load main configuration
    config = load_config(config_path)
    
    # Try to load advanced settings if available
    advanced_settings_path = joinpath(dirname(config_path), "advanced_settings.yaml")
    if isfile(advanced_settings_path)
        try
            advanced_settings = YAML.load_file(advanced_settings_path)
            config["advanced_settings"] = advanced_settings
            @info "Loaded advanced settings from $(advanced_settings_path)"
            
            # Verify the source hierarchy is properly configured
            verify_source_hierarchy(advanced_settings)
        catch e
            @warn "Failed to load advanced settings: $e"
        end
    end
    
    # Try to load plots settings if available
    plots_settings_path = joinpath(dirname(config_path), "plots.yaml")
    if isfile(plots_settings_path)
        try
            plots_settings = YAML.load_file(plots_settings_path)
            config["plots_settings"] = plots_settings
            @info "Loaded plot settings from $(plots_settings_path)"
        catch e
            @warn "Failed to load plot settings: $e"
        end
    end
    
    # Try to load species settings if available (trying both casing variants)
    for species_filename in ["species.yaml", "Species.yaml"]
        species_settings_path = joinpath(dirname(config_path), species_filename)
        if isfile(species_settings_path)
            try
                species_settings = YAML.load_file(species_settings_path)
                config["species_settings"] = species_settings
                @info "Loaded species settings from $(species_settings_path)"
                break  # Stop after successfully loading the first available file
            catch e
                @warn "Failed to load species settings from $(species_settings_path): $e"
            end
        end
    end
    
    # Initialize database connection
    db_path = config["general"]["database_path"]
    
    # Ensure database directory exists
    db_dir = dirname(db_path)
    isdir(db_dir) || mkpath(db_dir)
    
    # Connect to database
    try
        db = connect_database(db_path)
        @info "Successfully connected to database at $(db_path)"
    catch e
        @error "Failed to connect to database: $e"
        @warn "Will proceed without database connection"
        db = nothing
    end
    
    # Get database settings from configuration and apply them if database is connected
    if db !== nothing && haskey(config, "advanced_settings") && haskey(config["advanced_settings"], "database")
        db_settings = config["advanced_settings"]["database"]
        
        # Apply configuration to database
        if haskey(db_settings, "cache_size")
            try
                DuckDB.execute(db, "PRAGMA memory_limit='$(db_settings["cache_size"])MB'")
                @info "Set database cache size to $(db_settings["cache_size"])MB"
            catch e
                @warn "Failed to set database cache size: $e"
            end
        end
        
        # Other database settings could be applied here
    end
    
    # Initialize logging
    log_level = config["general"]["log_level"]
    init_logger(log_level)
    
    # Initialize hierarchical traversal system
    initialize_hierarchical_traversal()
    
    # Log initialization
    @info "JThermodynamicsData initialized" config_path=config_path database=db_path
    
    # Make sure to return both config and database connection
    return config, db
end

"""
    verify_source_hierarchy(config)

Verify that the source hierarchy is properly configured.

# Arguments
- `config`: The configuration dictionary containing the hierarchy
"""
function verify_source_hierarchy(config)
    # Default hierarchy if one is not found in the configuration
    default_hierarchy = [
        "theoretical",
        "group_contribution",
        "statistical_thermodynamics",
        "benson-group",
        "quantum_statistical",
        "gri-mech",
        "chemkin",
        "nasa-cea",
        "janaf",
        "thermoml",
        "tde",
        "nist-webbook",
        "burcat",
        "atct"
    ]

    # Check if hierarchy is defined in config
    if !haskey(config, "hierarchy") || !haskey(config["hierarchy"], "sources")
        @warn "No hierarchy section found in advanced settings, using default hierarchy"
        return default_hierarchy
    end
    
    # Get the hierarchy from config
    hierarchy = config["hierarchy"]["sources"]
    @info "Loaded source hierarchy from configuration: $(hierarchy)"
    
    # Verify that all expected sources are in the hierarchy
    expected_sources = [
        "theoretical",
        "group_contribution",
        "statistical_thermodynamics", 
        "benson-group",
        "quantum_statistical", 
        "gri-mech",
        "chemkin",
        "nasa-cea",
        "janaf",
        "thermoml",
        "tde",
        "nist-webbook",
        "burcat",
        "atct"
    ]
    
    # Check for missing sources
    missing_sources = setdiff(expected_sources, hierarchy)
    if !isempty(missing_sources)
        @warn "Some expected sources are missing from the hierarchy: $missing_sources"
    end
    
    # Check for additional sources that aren't expected
    additional_sources = setdiff(hierarchy, expected_sources)
    if !isempty(additional_sources)
        @info "Additional sources found in hierarchy: $additional_sources"
    end
    
    # Verify that the highest priority sources are correct
    if hierarchy[end] != "atct"
        @warn "The highest priority source should be 'atct', but found '$(hierarchy[end])'"
    end
    
    if hierarchy[end-1] != "burcat"
        @warn "The second highest priority source should be 'burcat', but found '$(hierarchy[end-1])'"
    end
    
    # Return the verified hierarchy
    return hierarchy
end

"""
    initialize_hierarchical_traversal()

Initialize the hierarchical traversal system for thermodynamic data.
This ensures all necessary components are loaded and ready.
"""
function initialize_hierarchical_traversal()
    @info "Initializing hierarchical traversal system"
    
    # Check that all parsers are available
    parser_modules = [
        "nasa", "chemkin", "janaf", "thermoml", "burcat", 
        "atct", "tde", "nist"
    ]
    
    for parser in parser_modules
        parser_path = joinpath(dirname(dirname(@__FILE__)), "parsers", "$(parser).jl")
        if isfile(parser_path)
            @info "Parser module found: $parser"
        else
            @warn "Parser module not found: $parser"
        end
    end
    
    # Check that json_storage utility is available
    json_storage_path = joinpath(dirname(dirname(@__FILE__)), "utils", "json_storage.jl")
    if isfile(json_storage_path)
        @info "JSON storage module found"
    else
        @error "JSON storage module not found! Hierarchical traversal may not work properly."
    end
    
    @info "Hierarchical traversal system initialized"
end

"""
    shutdown()

Clean up resources used by the package.
"""
function shutdown()
    # Close any open database connections
    close_all_connections()
    
    # Close loggers
    # (Add any other cleanup needed)
    
    @info "JThermodynamicsData shutdown complete"
end