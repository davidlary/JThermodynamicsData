"""
Hierarchical Thermodynamic Data Plotting System

This module provides a completely rewritten implementation of the plotting
system that properly implements the hierarchical data source selection,
ensuring that the most accurate source is always used for each species.

All settings are read from configuration files, with no hardcoded values.
"""

using Plots
using JSON
using YAML
using Dates
using Colors
using Printf
using Statistics
using DataFrames
using DuckDB

# Function to convert integers to superscript for logarithmic tick labels
function superscript(n::Int)
    if n == 0
        return "⁰"
    end
    
    superscripts = Dict(
        '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
        '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
        '-' => '⁻'
    )
    
    str = string(n)
    result = ""
    for c in str
        result *= superscripts[c]
    end
    
    return result
end

"""
    load_plot_config()

Load plot configuration from the YAML files.
Returns a merged configuration dictionary with all settings.
"""
function load_plot_config()
    # Find plot config file
    project_dir = dirname(dirname(dirname(dirname(@__FILE__))))
    plots_yaml_path = joinpath(project_dir, "config", "plots.yaml")
    settings_yaml_path = joinpath(project_dir, "config", "settings.yaml")
    advanced_yaml_path = joinpath(project_dir, "config", "advanced_settings.yaml")
    
    # Default configuration (fallback)
    default_config = Dict(
        "general" => Dict(
            "size" => [1200, 900],
            "dpi" => 300,
            "margin" => 10,
            "right_margin" => 20,
            "title_font_size" => 12,
            "temperature_points" => 100,
            "use_log_scale" => true
        ),
        "line_styles" => Dict(
            "best_source" => Dict(
                "line_width" => 2.5,
                "alpha" => 1.0,
                "style" => "solid",
                "color" => "blue"
            ),
            "other_sources" => Dict(
                "line_width" => 1.5,
                "alpha" => 0.7,
                "style" => "dash",
                "color_strategy" => "priority_based"
            )
        ),
        "uncertainty" => Dict(
            "ribbon_alpha" => 0.2
        ),
        "priority_colors" => Dict(
            "color_value_base" => 0.7,
            "color_value_step" => 0.05
        ),
        "output" => Dict(
            "file_format" => "png",
            "save_plots" => true,
            "plots_directory" => "plots",
            "documentation_directory" => "output/docs"
        ),
        "properties" => Dict(
            "cp" => Dict("label" => "Cp (J/mol/K)", "position" => 1),
            "h" => Dict("label" => "H (kJ/mol)", "position" => 2),
            "s" => Dict("label" => "S (J/mol/K)", "position" => 3),
            "g" => Dict("label" => "G (kJ/mol)", "position" => 4)
        )
    )
    
    config = deepcopy(default_config)
    
    # Load settings.yaml if available
    if isfile(settings_yaml_path)
        try
            settings = YAML.load_file(settings_yaml_path)
            # Merge relevant parts from settings
            if haskey(settings, "data_sources")
                config["data_sources"] = settings["data_sources"]
            end
            if haskey(settings, "visualization")
                config["visualization"] = settings["visualization"]
            end
            if haskey(settings, "theoretical_calculation")
                config["theoretical_calculation"] = settings["theoretical_calculation"]
            end
        catch e
            @warn "Failed to load settings.yaml: $e"
        end
    end
    
    # Load advanced_settings.yaml if available
    if isfile(advanced_yaml_path)
        try
            advanced = YAML.load_file(advanced_yaml_path)
            # Merge relevant parts
            if haskey(advanced, "hierarchy")
                config["hierarchy"] = advanced["hierarchy"]
            end
            if haskey(advanced, "uncertainty")
                config["advanced_uncertainty"] = advanced["uncertainty"]
            end
        catch e
            @warn "Failed to load advanced_settings.yaml: $e"
        end
    end
    
    # Load plots.yaml if available
    if isfile(plots_yaml_path)
        try
            plots_config = YAML.load_file(plots_yaml_path)
            # Merge all plot settings with priority to plots.yaml
            for (key, value) in plots_config
                config[key] = value
            end
        catch e
            @warn "Failed to load plots.yaml: $e"
        end
    end
    
    return config
end

"""
    get_temperature_range(config::Dict)

Get the temperature range from configuration.
"""
function get_temperature_range(config::Dict)
    # Try to get from species.yaml if in config
    if haskey(config, "species_settings") && 
       haskey(config["species_settings"], "temperature_range")
        return config["species_settings"]["temperature_range"]
    end
    
    # Try to get from general settings
    if haskey(config, "general") && 
       haskey(config["general"], "temperature_range")
        return config["general"]["temperature_range"]
    end
    
    # Default range
    return [100.0, 10000.0]
end

"""
    get_source_priority(config::Dict, source_name::String)

Get the priority for a data source.
"""
function get_source_priority(config::Dict, source_name::String)
    # For theoretical sources, check theoretical_calculation
    if startswith(lowercase(source_name), "theoretical") || 
       source_name in ["group_contribution", "statistical_thermodynamics", 
                      "quantum_statistical", "stat-thermo", "benson-group"]
        if haskey(config, "theoretical_calculation") && 
           haskey(config["theoretical_calculation"], "methods")
            for method in config["theoretical_calculation"]["methods"]
                if lowercase(method["name"]) == lowercase(source_name) || 
                   (startswith(lowercase(source_name), "theoretical") && 
                    lowercase(method["name"]) == "theoretical")
                    return get(method, "priority", 0)
                end
            end
        end
        # Default for theoretical methods
        return 0
    end
    
    # For experimental sources, check data_sources
    if haskey(config, "data_sources")
        for source in config["data_sources"]
            if lowercase(source["name"]) == lowercase(source_name)
                return get(source, "priority", 5)
            end
        end
    end
    
    # Default priority (middle range)
    return 5
end

"""
    get_source_reliability(config::Dict, source_name::String)

Get the reliability score for a data source.
"""
function get_source_reliability(config::Dict, source_name::String)
    # For theoretical sources
    if startswith(lowercase(source_name), "theoretical") || 
       source_name in ["group_contribution", "statistical_thermodynamics", 
                      "quantum_statistical", "stat-thermo", "benson-group"]
        # Default reliability for theoretical methods (0.0-3.0)
        return 2.5
    end
    
    # For experimental sources, check data_sources
    if haskey(config, "data_sources")
        for source in config["data_sources"]
            if lowercase(source["name"]) == lowercase(source_name)
                return get(source, "reliability_score", 4.0)
            end
        end
    end
    
    # Default reliability (middle range)
    return 4.0
end

"""
    standardize_source_name(source_name::String, priority::Int, config::Dict)

Standardize the display name for a data source.
"""
function standardize_source_name(source_name::String, priority::Integer, config::Dict)
    # Check if theoretical (priority <= 4 or name starts with theoretical)
    if priority <= 4 || startswith(lowercase(source_name), "theoretical") || 
       source_name in ["group_contribution", "statistical_thermodynamics", 
                      "quantum_statistical", "stat-thermo", "benson-group",
                      "machine_learning"]
        # Return specific theoretical source names for better differentiation
        if source_name == "theoretical"
            return "Theoretical"
        elseif source_name == "group_contribution" || source_name == "benson-group"
            return "Group-Contribution"
        elseif source_name == "statistical_thermodynamics" || source_name == "stat-thermo"
            return "Stat-Thermo"
        elseif source_name == "quantum_statistical"
            return "Quantum-Statistical"
        elseif source_name == "machine_learning"
            return "ML-Prediction"
        else
            return "Theoretical-" * source_name
        end
    end
    
    # Experimental sources - use proper capitalization
    words = split(source_name, "-")
    capitalized = [uppercase(first(word)) * lowercase(SubString(word, 2:lastindex(word))) 
                  for word in words]
    return join(capitalized, "-")
end

"""
    determine_source_color(source_name::String, priority::Int, is_best::Bool, config::Dict)

Determine the color for a data source based on configuration.
"""
function determine_source_color(source_name::String, priority::Integer, is_best::Bool, config::Dict)
    # Get color settings from config
    best_source_config = get(config, "line_styles", Dict())
    best_source_config = get(best_source_config, "best_source", Dict())
    best_color = get(best_source_config, "color", "blue")
    
    other_sources_config = get(config, "line_styles", Dict())
    other_sources_config = get(other_sources_config, "other_sources", Dict())
    color_strategy = get(other_sources_config, "color_strategy", "priority_based")
    
    # For best source, use configured color
    if is_best
        return Symbol(best_color)
    end
    
    # For other sources, use strategy from config
    if color_strategy == "priority_based"
        # Get color mapping settings
        priority_colors_config = get(config, "priority_colors", Dict())
        color_value_base = get(priority_colors_config, "color_value_base", 0.7)
        color_value_step = get(priority_colors_config, "color_value_step", 0.05)
        
        # Use distinct colors for different source types rather than just grayscale
        if priority <= 4  # Theoretical sources
            # Use different colors for different theoretical sources
            if source_name == "theoretical"
                return :gray50
            elseif source_name == "group_contribution" || source_name == "benson-group"
                return :brown
            elseif source_name == "statistical_thermodynamics" || source_name == "stat-thermo"
                return :purple
            elseif source_name == "quantum_statistical"
                return :teal
            elseif source_name == "machine_learning"
                return :green
            else
                return :gray30
            end
        else  # Experimental sources
            # Scale color by priority - higher priority = darker color
            priority_level = max(5, min(priority, 13)) - 5  # Normalize to 0-8 range
            color_value = color_value_base - (priority_level * color_value_step)
            return RGB(color_value, color_value, color_value)
        end
    elseif color_strategy == "source_based"
        # Fixed color per source with more vibrant, distinct colors
        source_colors = Dict(
            "ATCT" => :royalblue,
            "BURCAT" => :forestgreen,
            "NIST-WEBBOOK" => :crimson,
            "TDE" => :mediumorchid,
            "THERMOML" => :darkgoldenrod,
            "JANAF" => :darkorange,
            "NASA-CEA" => :deepskyblue,
            "CHEMKIN" => :hotpink,
            "GRI-MECH" => :gold,
            "Theoretical" => :slategray,
            "Group-Contribution" => :sienna,  
            "Stat-Thermo" => :mediumpurple,
            "Quantum-Statistical" => :teal,
            "ML-Prediction" => :limegreen
        )
        
        display_name = standardize_source_name(source_name, priority, config)
        return get(source_colors, display_name, :gray)
    else
        # Default grayscale
        return :gray
    end
end

"""
    determine_source_style(source_name::String, priority::Int, is_best::Bool, index::Int, config::Dict)

Determine the line style for a data source.
"""
function determine_source_style(source_name::String, priority::Integer, is_best::Bool, index::Integer, config::Dict)
    # Get style settings from config
    best_source_config = get(config, "line_styles", Dict())
    best_source_config = get(best_source_config, "best_source", Dict())
    best_style = Symbol(get(best_source_config, "style", "solid"))
    
    # For best source, always use solid line
    if is_best
        return best_style
    end
    
    # Alternate between different line styles
    styles = [:dash, :dot, :dashdot, :dashdotdot]
    style_index = (index - 1) % length(styles) + 1
    return styles[style_index]
end

"""
    get_hierarchical_sources(db::DuckDB.DB, species_name::String, config::Dict, temperature_range::Vector{Float64})

Get all available data sources for a species, properly sorted by hierarchical priority.
Returns a DataFrame with source information.
"""
function get_hierarchical_sources(db::DuckDB.DB, species_name::String, config::Dict, temperature_range::Vector{Float64})
    # Query all sources for this species within temperature range
    query = """
    SELECT 
        s.name AS species_name,
        s.formula,
        td.data_source,
        td.polynomial_type,
        td.temperature_min,
        td.temperature_max,
        td.reliability,
        td.data_json,
        td.uncertainty_json,
        ds.priority
    FROM 
        species s
    JOIN 
        thermodynamic_data td ON s.id = td.species_id
    JOIN 
        data_sources ds ON td.data_source = ds.name
    WHERE 
        s.name = ? AND
        td.temperature_min <= ? AND
        td.temperature_max >= ?
    GROUP BY 
        s.name, s.formula, td.data_source, td.polynomial_type, td.temperature_min,
        td.temperature_max, td.reliability, td.data_json, td.uncertainty_json, ds.priority
    ORDER BY 
        ds.priority ASC  -- Sort by priority ASCENDING (lowest to highest)
    """
    
    result = DuckDB.execute(db, query, [species_name, temperature_range[2], temperature_range[1]])
    sources_df = DataFrame(result)
    
    # If no sources found
    if nrow(sources_df) == 0
        @warn "No data sources found for $species_name in temperature range $temperature_range"
        return DataFrame()
    end
    
    # Check which sources are theoretical (priority <= 4 or name starts with "THEORETICAL")
    sources_df.is_theoretical = map(row -> 
        row.priority <= 4 || startswith(uppercase(row.data_source), "THEORETICAL"), 
        eachrow(sources_df))
    
    # Find all theoretical sources
    theoretical_sources = filter(row -> row.is_theoretical, sources_df)
    
    # Find all experimental sources
    experimental_sources = filter(row -> !row.is_theoretical, sources_df)
    
    # If we have theoretical sources, keep only the best one (highest priority)
    if nrow(theoretical_sources) > 0
        best_theoretical = sort(theoretical_sources, :priority, rev=true)[1:1, :]
    else
        best_theoretical = DataFrame()
    end
    
    # We need to handle two distinct scenarios:
    # 1. Only show the BEST theoretical approach (highest priority one)
    # 2. Show ALL actual experimental data sources
    
    # First, keep only the best theoretical source (highest priority)
    if nrow(theoretical_sources) > 0
        best_theoretical = sort(theoretical_sources, :priority, rev=true)[1:1, :]
    else
        best_theoretical = DataFrame()
    end
    
    # Now combine the best theoretical with ALL experimental sources
    combined_sources = vcat(best_theoretical, experimental_sources)
    
    # Sort by priority ascending to ensure proper hierarchy traversal
    sort!(combined_sources, :priority)
    
    return combined_sources
end

"""
    calculate_property_values(source_data::Dict, source_type::String, property::String, temperatures::Vector{Float64})

Calculate thermodynamic property values for a list of temperatures based on source data.
"""
function calculate_property_values(source_data::Dict, source_type::String, property::String, temperatures::Vector{Float64})
    # Parse the data JSON if needed
    data = source_data
    if typeof(source_data) == String
        try
            data = JSON.parse(source_data)
        catch e
            @warn "Failed to parse data JSON: $e"
            return Float64[], Int[]
        end
    end
    
    # Initialize results
    values = Float64[]
    valid_indices = Int[]
    
    # Handle different polynomial types
    if source_type == "nasa7" || source_type == "nasa9"
        # Extract coefficients and temperature ranges based on actual data structure
        # Check if the data has the expected structure
        if haskey(data, "low_temp") && haskey(data, "high_temp")
            # Structure with low_temp and high_temp keys
            coeffs = [data["low_temp"]["coefficients"], data["high_temp"]["coefficients"]]
            temp_ranges = [data["low_temp"]["range"], data["high_temp"]["range"]]
        elseif haskey(data, "coefficients") && haskey(data, "temperature_ranges")
            # Alternative structure with direct coefficients and temperature_ranges keys
            coeffs = data["coefficients"]
            temp_ranges = data["temperature_ranges"]
        else
            @warn "Unexpected data structure: missing required keys for NASA polynomials"
            return Float64[], Int[]
        end
        
        # Calculate property for each temperature
        for (i, T) in enumerate(temperatures)
            # Find applicable temperature range
            range_idx = 0
            for (idx, range) in enumerate(temp_ranges)
                if range[1] <= T && T <= range[2]
                    range_idx = idx
                    break
                end
            end
            
            # Skip if no applicable range
            if range_idx == 0
                continue
            end
            
            # Get coefficients for this range
            coefs = coeffs[range_idx]
            
            # Calculate property
            if property == "cp"
                # NASA-7 polynomial: Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
                cp_over_r = coefs[1] + coefs[2]*T + coefs[3]*T^2 + coefs[4]*T^3 + coefs[5]*T^4
                push!(values, cp_over_r * 8.31446)  # Convert to J/mol/K
                push!(valid_indices, i)
                
            elseif property == "h"
                # NASA-7 polynomial: H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
                h_over_rt = coefs[1] + coefs[2]*T/2 + coefs[3]*T^2/3 + coefs[4]*T^3/4 + coefs[5]*T^4/5 + coefs[6]/T
                push!(values, h_over_rt * 8.31446 * T / 1000.0)  # Convert to kJ/mol
                push!(valid_indices, i)
                
            elseif property == "s"
                # NASA-7 polynomial: S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
                s_over_r = coefs[1]*log(T) + coefs[2]*T + coefs[3]*T^2/2 + coefs[4]*T^3/3 + coefs[5]*T^4/4 + coefs[7]
                push!(values, s_over_r * 8.31446)  # Convert to J/mol/K
                push!(valid_indices, i)
                
            elseif property == "g"
                # Calculate G = H - T*S
                h_over_rt = coefs[1] + coefs[2]*T/2 + coefs[3]*T^2/3 + coefs[4]*T^3/4 + coefs[5]*T^4/5 + coefs[6]/T
                s_over_r = coefs[1]*log(T) + coefs[2]*T + coefs[3]*T^2/2 + coefs[4]*T^3/3 + coefs[5]*T^4/4 + coefs[7]
                
                h = h_over_rt * 8.31446 * T / 1000.0  # kJ/mol
                s = s_over_r * 8.31446 / 1000.0       # kJ/mol/K
                g = h - T * s
                
                push!(values, g)  # kJ/mol
                push!(valid_indices, i)
            end
        end
        
    elseif source_type == "tabular"
        # For tabular data, we need to interpolate
        if !haskey(data, "data_points")
            return Float64[], Int[]
        end
        
        data_points = data["data_points"]
        
        # Sort data points by temperature
        sort!(data_points, by = p -> p["temperature"])
        
        # Extract temperatures and property values for interpolation
        table_temps = [p["temperature"] for p in data_points]
        
        # Property mapping
        prop_map = Dict(
            "cp" => "Cp", 
            "h" => "H", 
            "s" => "S", 
            "g" => "G"
        )
        
        # Check if property exists
        prop_key = prop_map[property]
        if !all(haskey(p, prop_key) for p in data_points)
            return Float64[], Int[]
        end
        
        table_values = [p[prop_key] for p in data_points]
        
        # Calculate for each temperature using linear interpolation
        for (i, T) in enumerate(temperatures)
            # Skip if outside range
            if T < table_temps[1] || T > table_temps[end]
                continue
            end
            
            # Find interpolation points
            idx_low = findlast(t -> t <= T, table_temps)
            idx_high = findfirst(t -> t >= T, table_temps)
            
            # Exact match
            if table_temps[idx_low] == T
                push!(values, table_values[idx_low])
                push!(valid_indices, i)
            # Interpolate
            elseif idx_low != idx_high
                t_low = table_temps[idx_low]
                t_high = table_temps[idx_high]
                v_low = table_values[idx_low]
                v_high = table_values[idx_high]
                
                # Linear interpolation
                factor = (T - t_low) / (t_high - t_low)
                interp_value = v_low + factor * (v_high - v_low)
                
                push!(values, interp_value)
                push!(valid_indices, i)
            end
        end
    end
    
    return values, valid_indices
end

"""
    plot_species_with_hierarchy(species_name::String, db::DuckDB.DB, config::Dict; save_plot::Bool=true, display_plot::Bool=true)

Create a comprehensive thermodynamic plot for a species, showing all available sources
with proper hierarchical selection and uncertainty visualization.
"""
function plot_species_with_hierarchy(species_name::String, db::DuckDB.DB, config::Dict; 
                                    save_plot::Bool=true, display_plot::Bool=true)
    # Get global settings from config
    temp_range = get_temperature_range(config)
    general_settings = get(config, "general", Dict())
    size_settings = get(general_settings, "size", [1200, 900])
    dpi_settings = get(general_settings, "dpi", 300)
    margin_settings = get(general_settings, "margin", 10)
    right_margin_settings = get(general_settings, "right_margin", 20)
    use_log_scale = get(general_settings, "use_log_scale", true)
    temperature_points = get(general_settings, "temperature_points", 100)
    
    # Generate temperature points
    if use_log_scale
        # Logarithmic temperature scale
        temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=temperature_points)
    else
        # Linear temperature scale
        temps = range(temp_range[1], temp_range[2], length=temperature_points)
    end
    
    # Get all sources for this species with hierarchical order
    # Convert temperature range to Float64 explicitly to avoid type issues
    float_temp_range = Float64[Float64(temp_range[1]), Float64(temp_range[2])]
    sources_df = get_hierarchical_sources(db, species_name, config, float_temp_range)
    
    if nrow(sources_df) == 0
        @warn "No data sources available for $species_name in temperature range $temp_range"
        return nothing
    end
    
    # Extract formula if available
    formula = sources_df.formula[1]
    
    # Determine best source (highest priority)
    best_source_row = sources_df[end, :]  # Last row has highest priority
    best_source_name = best_source_row.data_source
    best_source_priority = best_source_row.priority
    
    @info "Best source for $species_name: $best_source_name (priority: $best_source_priority)"
    
    # Get display names for all sources
    source_display_names = Dict()
    for row in eachrow(sources_df)
        source_display_names[row.data_source] = standardize_source_name(
            row.data_source, row.priority, config
        )
    end
    
    # Create property settings
    properties_config = get(config, "properties", Dict())
    properties = [
        (get(get(properties_config, "cp", Dict()), "position", 1), "cp", get(get(properties_config, "cp", Dict()), "label", "Cp (J/mol/K)")),
        (get(get(properties_config, "h", Dict()), "position", 2), "h", get(get(properties_config, "h", Dict()), "label", "H (kJ/mol)")),
        (get(get(properties_config, "s", Dict()), "position", 3), "s", get(get(properties_config, "s", Dict()), "label", "S (J/mol/K)")),
        (get(get(properties_config, "g", Dict()), "position", 4), "g", get(get(properties_config, "g", Dict()), "label", "G (kJ/mol)"))
    ]
    
    # Create the plot with 2x2 layout
    plt = plot(
        layout=(2,2), 
        size=(size_settings[1], size_settings[2]), 
        dpi=dpi_settings, 
        margin=margin_settings * Plots.mm, 
        right_margin=right_margin_settings * Plots.mm,
        legend=:outerright
    )
    
    # Track sources we've already plotted (avoid duplicates)
    processed_sources = Dict()
    
    # Get uncertainty settings
    uncertainty_settings = get(config, "uncertainty", Dict())
    ribbon_alpha = get(uncertainty_settings, "ribbon_alpha", 0.2)
    
    # Process each source in hierarchical order (lowest to highest priority)
    for (idx, row) in enumerate(eachrow(sources_df))
        source_name = row.data_source
        
        # Skip if already processed
        if haskey(processed_sources, source_name)
            continue
        end
        processed_sources[source_name] = true
        
        @info "  Processing source: $source_name (priority: $(row.priority))"
        
        # Determine if this is the best source
        is_best = source_name == best_source_name
        
        # Get display name
        display_name = source_display_names[source_name]
        
        # Set line style properties
        line_color = determine_source_color(source_name, row.priority, is_best, config)
        line_style = determine_source_style(source_name, row.priority, is_best, idx, config)
        
        # Set line width and alpha based on config
        line_styles = get(config, "line_styles", Dict())
        if is_best
            best_styles = get(line_styles, "best_source", Dict())
            line_width = get(best_styles, "line_width", 2.5)
            alpha = get(best_styles, "alpha", 1.0)
        else
            other_styles = get(line_styles, "other_sources", Dict())
            line_width = get(other_styles, "line_width", 1.5)
            alpha = get(other_styles, "alpha", 0.7)
        end
        
        # Parse uncertainty if available
        uncertainty_value = 0.0
        if row.uncertainty_json !== missing && !isempty(row.uncertainty_json)
            try
                uncertainty_data = JSON.parse(row.uncertainty_json)
                uncertainty_value = get(uncertainty_data, "value", 0.1)
            catch e
                @warn "Failed to parse uncertainty JSON for $source_name: $e"
                uncertainty_value = 0.1  # Default uncertainty (10%)
            end
        end
        
        # Calculate and plot properties
        for (idx, property, ylabel) in properties
            # Calculate property values and uncertainties
            values = Float64[]
            valid_indices = Int[]
            
            # Handle string vs dict data
            data_obj = row.data_json
            if typeof(data_obj) == String
                try
                    data_obj = JSON.parse(data_obj)
                catch e
                    @warn "Failed to parse data JSON: $e"
                    continue
                end
            end
            
            values, valid_indices = calculate_property_values(
                data_obj, row.polynomial_type, property, temps
            )
            
            # Skip if no valid values
            if isempty(values)
                @warn "  No valid values for $property from $source_name"
                continue
            end
            
            # Apply uncertainty to values
            uncertainties = values .* uncertainty_value
            
            # Filter temperatures for valid indices
            valid_temps = temps[valid_indices]
            
            # Generate logarithmic tick marks if needed
            if use_log_scale
                # Create logarithmic ticks for each decade
                log_min = log10(temp_range[1])
                log_max = log10(temp_range[2])
                decades = floor(Int, log_min):ceil(Int, log_max)
                
                # Generate ticks within each decade
                tick_positions = Float64[]
                tick_labels = String[]
                
                for decade in decades
                    # Major tick for the decade
                    push!(tick_positions, 10.0^decade)
                    push!(tick_labels, "10$(superscript(decade))")
                    
                    # Minor ticks within the decade (2-9)
                    for i in 2:9
                        minor_tick = i * 10.0^decade
                        if minor_tick >= temp_range[1] && minor_tick <= temp_range[2]
                            push!(tick_positions, minor_tick)
                            # Use empty string for minor ticks
                            push!(tick_labels, "")
                        end
                    end
                end
                
                # Plot for this property
                if is_best
                    # Best source - highlight with thicker line and uncertainty ribbon
                    plot!(
                        plt, 
                        valid_temps, 
                        values, 
                        subplot=idx,
                        label=display_name,
                        xlabel="Temperature (K)",
                        ylabel=ylabel,
                        title="$(uppercase(property)) for $species_name",
                        xscale=:log10,
                        xticks=(tick_positions, tick_labels),
                        xminorticks=10,
                        line=(:blue, 4.0, :solid),  # Force extra thick blue line
                        ribbon=(values .- (values .- uncertainties), (values .+ uncertainties) .- values),
                        fillalpha=0.3,
                        fillcolor=:lightblue
                    )
                else
                    # Other sources - no uncertainty ribbon
                    plot!(
                        plt, 
                        valid_temps, 
                        values, 
                        subplot=idx,
                        label=display_name,
                        line=(line_color, line_width, line_style),
                        alpha=alpha,
                        xscale=:log10
                    )
                end
            else
                # Linear scale plotting
                if is_best
                    # Best source
                    plot!(
                        plt, 
                        valid_temps, 
                        values, 
                        subplot=idx,
                        label=display_name,
                        xlabel="Temperature (K)",
                        ylabel=ylabel,
                        title="$(uppercase(property)) for $species_name",
                        line=(:blue, 4.0, :solid),
                        ribbon=(values .- (values .- uncertainties), (values .+ uncertainties) .- values),
                        fillalpha=0.3,
                        fillcolor=:lightblue
                    )
                else
                    # Other sources - no uncertainty ribbon
                    plot!(
                        plt, 
                        valid_temps, 
                        values, 
                        subplot=idx,
                        label=display_name,
                        line=(line_color, line_width, line_style),
                        alpha=alpha
                    )
                end
            end
        end
    end
    
    # Add a main title
    title_font_size = get(general_settings, "title_font_size", 12)
    plot!(plt, title="Thermodynamic Properties for $species_name ($formula)", 
          titleloc=:center, titlefontsize=title_font_size)
    
    # Save the plot if requested
    if save_plot
        output_settings = get(config, "output", Dict())
        plots_dir = get(output_settings, "plots_directory", "plots")
        file_format = get(output_settings, "file_format", "png")
        
        # Ensure directory exists
        mkpath(plots_dir)
        
        savefig(plt, joinpath(plots_dir, "$(species_name)_all_sources.$(file_format)"))
        @info "Plot saved to $(plots_dir)/$(species_name)_all_sources.$(file_format)"
        
        # Create documentation
        docs_dir = get(output_settings, "documentation_directory", "output/docs")
        mkpath(docs_dir)
        
        # Create document with source information
        # Convert temperature range to Float64 explicitly to match function signature
        float_temp_range = Float64[Float64(temp_range[1]), Float64(temp_range[2])]
        create_source_documentation(species_name, sources_df, best_source_name, float_temp_range, 
                                   config, joinpath(docs_dir, "$(species_name)_sources.md"))
    end
    
    # Return the plot
    return plt
end

"""
    create_source_documentation(species_name::String, sources_df::DataFrame, best_source::String, 
                               temp_range::Vector{Float64}, config::Dict, output_path::String)

Create Markdown documentation about the data sources used for a species.
"""
function create_source_documentation(species_name::String, sources_df::DataFrame, 
                                    best_source::String, temp_range::Vector{Float64}, 
                                    config::Dict, output_path::String)
    # Extract formula if available
    formula = sources_df.formula[1]
    
    # Get best source row
    best_source_row = filter(row -> row.data_source == best_source, sources_df)[1, :]
    
    # Create standardized display names
    source_display_names = Dict()
    for row in eachrow(sources_df)
        source_display_names[row.data_source] = standardize_source_name(
            row.data_source, row.priority, config
        )
    end
    
    # Create documentation
    open(output_path, "w") do f
        write(f, "# Thermodynamic Data Source for $species_name ($formula)\n\n")
        
        # Selected source section
        write(f, "## Selected Source\n")
        write(f, "- **Source:** $(source_display_names[best_source])\n")
        write(f, "- **Priority:** $(best_source_row.priority)\n")
        write(f, "- **Reliability score:** $(best_source_row.reliability)\n")
        
        # Parse uncertainty if available
        if best_source_row.uncertainty_json !== missing && !isempty(best_source_row.uncertainty_json)
            try
                uncertainty_data = JSON.parse(best_source_row.uncertainty_json)
                uncertainty_value = get(uncertainty_data, "value", 0.1)
                write(f, "- **Uncertainty:** $(uncertainty_value * 100)%\n")
            catch
                write(f, "- **Uncertainty:** Unknown\n")
            end
        else
            write(f, "- **Uncertainty:** Unknown\n")
        end
        
        write(f, "- **Temperature range:** $(temp_range[1]) - $(temp_range[2]) K\n")
        write(f, "- **Polynomial type:** $(best_source_row.polynomial_type)\n\n")
        
        # Add information about all available sources
        write(f, "\n## All Available Sources\n")
        write(f, "| Source | Priority | Reliability | Temperature Range |\n")
        write(f, "|--------|----------|-------------|-------------------|\n")
        
        # Sort sources by priority for the table (highest to lowest)
        sorted_sources = sort(sources_df, :priority, rev=true)
        
        for row in eachrow(sorted_sources)
            display_name = source_display_names[row.data_source]
            temp_range_str = "$(Int(row.temperature_min))-$(Int(row.temperature_max)) K"
            write(f, "| $(display_name) | $(row.priority) | $(row.reliability) | $(temp_range_str) |\n")
        end
        
        # Add hierarchy information
        write(f, "\n## Hierarchical Selection Process\n")
        write(f, "Sources were evaluated in hierarchical order from lowest to highest priority:\n\n")
        
        for (i, row) in enumerate(eachrow(sources_df))
            display_name = source_display_names[row.data_source]
            write(f, "$(i). $(display_name) (Priority: $(row.priority))\n")
        end
        
        write(f, "\nFinal selected source: **$(source_display_names[best_source])** (highest priority available)\n")
        
        # Add generation information
        write(f, "\n---\n")
        write(f, "*Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) using JThermodynamicsData hierarchical selection system*\n")
    end
    
    return output_path
end

"""
    run_comprehensive_plotting(config_path::String=""; species_list::Vector{String}=String[])

Run the comprehensive plotting for all species or a specific list.
"""
function run_comprehensive_plotting(config_path::String=""; species_list::Vector{String}=String[])
    # Load configuration
    if isempty(config_path)
        config_path = joinpath(dirname(dirname(dirname(dirname(@__FILE__)))), "config", "settings.yaml")
    end
    
    # Load config and connect to database
    config, db = initialize(config_path)
    if db === nothing
        error("Failed to connect to database")
    end
    
    # Load plot configuration
    plot_config = load_plot_config()
    
    # Merge configs
    merged_config = deepcopy(config)
    for (k, v) in plot_config
        merged_config[k] = v
    end
    
    # If no species list provided, get from configuration
    if isempty(species_list)
        if haskey(merged_config, "species_settings") && haskey(merged_config["species_settings"], "species")
            species_list = merged_config["species_settings"]["species"]
        else
            # Query database for all species
            result = DuckDB.execute(db, "SELECT DISTINCT name FROM species ORDER BY name")
            df = DataFrame(result)
            species_list = df.name
        end
    end
    
    # Start timer
    start_time = now()
    
    # Process each species
    successful_species = String[]
    for (i, species) in enumerate(species_list)
        @info "Processing $(i)/$(length(species_list)): $(species)"
        try
            plt = plot_species_with_hierarchy(species, db, merged_config; save_plot=true, display_plot=false)
            if plt !== nothing
                push!(successful_species, species)
            end
        catch e
            @warn "Failed to process $species: $e"
        end
    end
    
    # Generate summary report
    create_summary_report(species_list, successful_species, db, merged_config)
    
    # Calculate elapsed time
    elapsed = now() - start_time
    
    @info "Completed processing $(length(successful_species))/$(length(species_list)) species in $(Dates.format(elapsed, "HH:MM:SS"))"
    
    return successful_species
end

"""
    create_summary_report(all_species::Vector{String}, success_species::Vector{String}, 
                         db::DuckDB.DB, config::Dict)

Create a summary report of all processed species.
"""
function create_summary_report(all_species::Vector{String}, success_species::Vector{String}, 
                              db::DuckDB.DB, config::Dict)
    @info "Creating summary report..."
    
    # Output directory
    output_settings = get(config, "output", Dict())
    output_dir = get(output_settings, "documentation_directory", "output/docs")
    output_dir = dirname(output_dir)
    mkpath(output_dir)
    
    # Create source statistics
    source_counts = Dict{String, Int}()
    source_species = Dict{String, Vector{String}}()
    
    for species in success_species
        # Get all sources for this species
        query = """
        SELECT 
            td.data_source,
            ds.priority,
            COUNT(*) as count
        FROM 
            species s
        JOIN 
            thermodynamic_data td ON s.id = td.species_id
        JOIN 
            data_sources ds ON td.data_source = ds.name
        WHERE 
            s.name = ?
        GROUP BY 
            td.data_source, ds.priority
        """
        
        result = DuckDB.execute(db, query, [species])
        sources_df = DataFrame(result)
        
        for row in eachrow(sources_df)
            source_name = row.data_source
            
            # Standardize display name
            display_name = standardize_source_name(source_name, row.priority, config)
            
            # Update counts
            source_counts[display_name] = get(source_counts, display_name, 0) + 1
            
            # Track species for each source
            if !haskey(source_species, display_name)
                source_species[display_name] = String[]
            end
            push!(source_species[display_name], species)
        end
    end
    
    # Create main summary report
    summary_file = joinpath(output_dir, "data_source_summary.md")
    open(summary_file, "w") do f
        write(f, "# Thermodynamic Data Source Summary\n\n")
        
        write(f, "## Overview\n")
        write(f, "- Total species: $(length(all_species))\n")
        write(f, "- Successfully processed: $(length(success_species))\n")
        write(f, "- Success rate: $(round(length(success_species) / length(all_species) * 100, digits=1))%\n\n")
        
        write(f, "## Source Statistics\n")
        write(f, "| Source | Count | Percentage |\n")
        write(f, "|--------|-------|------------|\n")
        
        total_sources = sum(values(source_counts))
        
        for (source, count) in sort(collect(source_counts), by=x->x[2], rev=true)
            percentage = round(count / length(success_species) * 100, digits=1)
            write(f, "| $(source) | $(count) | $(percentage)% |\n")
        end
        
        write(f, "\n## Best Source Distribution\n")
        write(f, "| Source | Species Count | Species List |\n")
        write(f, "|--------|---------------|-------------|\n")
        
        # Count best sources
        best_source_counts = Dict{String, Int}()
        best_source_species = Dict{String, Vector{String}}()
        
        for species in success_species
            # Get best source for this species
            query = """
            SELECT 
                s.name,
                td.data_source,
                ds.priority
            FROM 
                species s
            JOIN 
                thermodynamic_data td ON s.id = td.species_id
            JOIN 
                data_sources ds ON td.data_source = ds.name
            WHERE 
                s.name = ?
            ORDER BY 
                ds.priority DESC
            LIMIT 1
            """
            
            result = DuckDB.execute(db, query, [species])
            df = DataFrame(result)
            
            if nrow(df) > 0
                source_name = df[1, :data_source]
                priority = df[1, :priority]
                
                # Standardize display name
                display_name = standardize_source_name(source_name, priority, config)
                
                # Update counts
                best_source_counts[display_name] = get(best_source_counts, display_name, 0) + 1
                
                # Track species for each source
                if !haskey(best_source_species, display_name)
                    best_source_species[display_name] = String[]
                end
                push!(best_source_species[display_name], species)
            end
        end
        
        # Write best source statistics
        for (source, count) in sort(collect(best_source_counts), by=x->x[2], rev=true)
            species_list = join(sort(best_source_species[source])[1:min(5, length(best_source_species[source]))], ", ")
            if length(best_source_species[source]) > 5
                species_list *= ", ... ($(length(best_source_species[source])) total)"
            end
            write(f, "| $(source) | $(count) | $(species_list) |\n")
        end
        
        write(f, "\n## Species List\n")
        write(f, "| Species | Processed | Sources |\n")
        write(f, "|---------|-----------|--------|\n")
        
        for species in sort(all_species)
            processed = species in success_species ? "✓" : "❌"
            
            # Count sources if processed
            sources_count = 0
            if species in success_species
                query = """
                SELECT COUNT(DISTINCT td.data_source) as count
                FROM species s
                JOIN thermodynamic_data td ON s.id = td.species_id
                WHERE s.name = ?
                """
                
                result = DuckDB.execute(db, query, [species])
                df = DataFrame(result)
                
                if nrow(df) > 0
                    sources_count = df[1, :count]
                end
            end
            
            write(f, "| $(species) | $(processed) | $(sources_count) |\n")
        end
        
        write(f, "\n---\n")
        write(f, "*Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) using JThermodynamicsData hierarchical selection system*\n")
    end
    
    @info "Summary report saved to $summary_file"
    return summary_file
end

# Export functions
export plot_species_with_hierarchy, run_comprehensive_plotting, load_plot_config