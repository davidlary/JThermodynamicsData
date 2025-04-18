module JThermodynamicsData

using YAML
using DuckDB
using DataFrames
using CSV
using HTTP
using JSON
using Measurements
using MonteCarloMeasurements
using Interpolations
using Statistics
using Plots
using Unitful
using UnitfulMoles
using LightXML
using StaticArrays
using SpecialFunctions

# Core functionality
include("JThermodynamicsData/core/config.jl")
include("JThermodynamicsData/core/types.jl")
include("JThermodynamicsData/core/constants.jl")

# Database functionality
include("JThermodynamicsData/database/schema.jl")
include("JThermodynamicsData/database/connection.jl")
include("JThermodynamicsData/database/queries.jl")

# Data parsers for different formats
include("JThermodynamicsData/parsers/nasa.jl")
include("JThermodynamicsData/parsers/chemkin.jl")
include("JThermodynamicsData/parsers/janaf.jl")
include("JThermodynamicsData/parsers/thermoml.jl")
include("JThermodynamicsData/parsers/burcat.jl")
include("JThermodynamicsData/parsers/atct.jl")
include("JThermodynamicsData/parsers/tde.jl")

# Theoretical models
include("JThermodynamicsData/models/group_contribution.jl")
include("JThermodynamicsData/models/statistical_thermodynamics.jl")
include("JThermodynamicsData/models/machine_learning.jl")
include("JThermodynamicsData/models/quantum_chemistry.jl")

# Utility functions
include("JThermodynamicsData/utils/fetcher.jl")
include("JThermodynamicsData/utils/cache.jl")
include("JThermodynamicsData/utils/converter.jl")
include("JThermodynamicsData/utils/uncertainty.jl")
include("JThermodynamicsData/utils/logger.jl")
include("JThermodynamicsData/utils/refinement.jl")

# API and interface
include("JThermodynamicsData/api/query.jl")
include("JThermodynamicsData/api/rest.jl")

# Visualization
include("JThermodynamicsData/visualization/plots.jl")
include("JThermodynamicsData/visualization/comparison.jl")

# Export main functions
export initialize, query_species, query_properties
export calculate_properties, plot_properties
export update_database, get_uncertainty
export list_available_species, list_available_databases
export ThermodynamicData, ThermodynamicProperty

# Export configuration functions
export load_config, validate_config, merge_configs, save_config

# Export database functions
export init_database, close_database

# Export parser functions 
export nasa_to_database, nasa9_to_database, burcat_to_database, chemkin_to_database
export janaf_to_database, thermoml_to_database, tde_to_database, atct_to_database
export parse_nasa7_file, parse_nasa9_file, parse_burcat_file, parse_chemkin_file
export parse_janaf_file, parse_thermoml_file, parse_tde_file, parse_atct_file

# Export thermodynamic models
export estimate_properties_group_contribution, estimate_properties_statistical_thermodynamics
export estimate_properties_quantum_chemistry, estimate_properties_machine_learning
export benson_group_additivity, joback_method

# Export refinement functions
export progressively_refine_thermodynamic_data, get_all_thermodynamic_data_for_species
export generate_source_usage_report, create_markdown_documentation, create_json_documentation

# Export documentation functions
export create_markdown_documentation, create_json_documentation

# Export utility functions
export fetch_data_source

# Main initialization function
function __init__()
    # Load configuration if config file exists
    config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
    if isfile(config_path)
        global_config = load_config(config_path)
    end
end

# Main entry point for programmatic usage
function initialize(config_path::String=joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml"))
    # Load configuration
    config = load_config(config_path)
    
    # Initialize database
    db = init_database(config)
    
    # Create cache directory if it doesn't exist
    cache_dir = config["general"]["cache_directory"]
    mkpath(cache_dir)
    
    # Initialize logger
    init_logger(config["general"]["log_level"])
    
    # Return initialized context
    return config, db
end

end # module