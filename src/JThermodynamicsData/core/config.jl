"""
    load_config(config_path::String)

Load configuration from a YAML file.
"""
function load_config(config_path::String)
    if !isfile(config_path)
        error("Configuration file not found: $config_path")
    end
    
    config = YAML.load_file(config_path)
    validate_config(config)
    
    return config
end

"""
    validate_config(config::Dict)

Validate the configuration settings.
"""
function validate_config(config::Dict)
    # Required top-level keys
    required_keys = ["general", "data_sources", "database", "theoretical_calculation", "uncertainty"]
    
    for key in required_keys
        if !haskey(config, key)
            error("Missing required configuration section: $key")
        end
    end
    
    # Validate temperature range
    temp_range = config["general"]["temperature_range"]
    if !(temp_range isa Vector) || length(temp_range) != 2 || 
       temp_range[1] >= temp_range[2] || temp_range[1] <= 0
        error("Invalid temperature range: $temp_range. Must be [min, max] with min > 0 and max > min.")
    end
    
    # Validate data sources
    for (i, source) in enumerate(config["data_sources"])
        required_source_keys = ["name", "format", "enabled", "priority"]
        for key in required_source_keys
            if !haskey(source, key)
                error("Missing required key '$key' in data source #$i: $(get(source, "name", "unnamed"))")
            end
        end
        
        # Check if either url or path is provided
        if !haskey(source, "url") && !haskey(source, "path")
            error("Either 'url' or 'path' must be provided for data source: $(source["name"])")
        end
    end
    
    # Sort data sources by priority
    sort!(config["data_sources"], by = s -> s["priority"])
    
    return true
end

"""
    merge_configs(base_config::Dict, override_config::Dict)

Merge two configurations, with override_config taking precedence.
"""
function merge_configs(base_config::Dict, override_config::Dict)
    merged = deepcopy(base_config)
    
    for (key, value) in override_config
        if value isa Dict && haskey(merged, key) && merged[key] isa Dict
            merged[key] = merge_configs(merged[key], value)
        else
            merged[key] = value
        end
    end
    
    return merged
end

"""
    save_config(config::Dict, config_path::String)

Save configuration to a YAML file.
"""
function save_config(config::Dict, config_path::String)
    YAML.write_file(config_path, config)
    return config_path
end