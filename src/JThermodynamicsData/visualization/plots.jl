"""
Functions for plotting thermodynamic data.
"""

"""
    plot_property(property_data::Dict, property_name::String, config::Dict)

Plot a thermodynamic property as a function of temperature.
"""
function plot_property(property_data::Dict, property_name::String, config::Dict)
    # Extract data
    temperatures = property_data["temperatures"]
    values = property_data["properties"][property_name]["values"]
    uncertainties = property_data["properties"][property_name]["uncertainties"]
    units = property_data["properties"][property_name]["units"]
    species_name = property_data["species_name"]
    
    # Set up plot theme
    theme = config["visualization"]["default_theme"]
    default_theme = theme == "default" ? :default : Symbol(theme)
    theme!(default_theme)
    
    # Create plot
    p = plot(
        temperatures, 
        values, 
        xlabel="Temperature (K)", 
        ylabel="$property_name ($(units))",
        title="$property_name for $species_name",
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Add uncertainty visualization if enabled
    if config["visualization"]["show_uncertainty"]
        uncertainty_style = config["visualization"]["uncertainty_style"]
        
        if uncertainty_style == "ribbon"
            # Ribbon style (shaded area)
            plot!(
                p,
                temperatures,
                values,
                ribbon=uncertainties,
                fillalpha=0.3
            )
        elseif uncertainty_style == "errorbar"
            # Error bars
            scatter!(
                p,
                temperatures[1:max(1, div(length(temperatures), 20)):end],  # Downsample for clarity
                values[1:max(1, div(length(values), 20)):end],
                yerror=uncertainties[1:max(1, div(length(uncertainties), 20)):end],
                markersize=3,
                markerstrokewidth=0.5
            )
        elseif uncertainty_style == "band"
            # Uncertainty bands (upper and lower curves)
            plot!(
                p,
                temperatures,
                values + uncertainties,
                linestyle=:dash,
                linecolor=:black,
                linewidth=1,
                alpha=0.5
            )
            
            plot!(
                p,
                temperatures,
                values - uncertainties,
                linestyle=:dash,
                linecolor=:black,
                linewidth=1,
                alpha=0.5
            )
        end
    end
    
    return p
end

"""
    plot_properties(property_data::Dict, config::Dict)

Plot all thermodynamic properties as a function of temperature.
"""
function plot_properties(property_data::Dict, config::Dict)
    # Create individual plots
    p_cp = plot_property(property_data, "Cp", config)
    p_h = plot_property(property_data, "H", config)
    p_s = plot_property(property_data, "S", config)
    p_g = plot_property(property_data, "G", config)
    
    # Combine into a 2x2 layout
    p = plot(p_cp, p_h, p_s, p_g, layout=(2, 2), size=(800, 600))
    
    return p
end

"""
    plot_species_comparison(conn::DuckDB.DB, species_names::Vector{String}, 
                          property_name::String, temperature_range::Vector{Float64}, 
                          config::Dict; step::Float64=10.0)

Compare a property across multiple species.
"""
function plot_species_comparison(conn::DuckDB.DB, species_names::Vector{String}, 
                               property_name::String, temperature_range::Vector{Float64}, 
                               config::Dict; step::Float64=10.0)
    # Set up plot theme
    theme = config["visualization"]["default_theme"]
    default_theme = theme == "default" ? :default : Symbol(theme)
    theme!(default_theme)
    
    # Create plot
    p = plot(
        xlabel="Temperature (K)", 
        ylabel="$property_name",
        title="$property_name Comparison",
        grid=true,
        framestyle=:box
    )
    
    # Calculate properties for each species
    for species_name in species_names
        try
            result = calculate_properties(conn, species_name, temperature_range, config, step=step)
            
            temperatures = result["temperatures"]
            values = result["properties"][property_name]["values"]
            uncertainties = result["properties"][property_name]["uncertainties"]
            units = result["properties"][property_name]["units"]
            
            # Add to plot
            plot!(
                p,
                temperatures, 
                values,
                label=species_name,
                lw=2
            )
            
            # Add uncertainty if enabled
            if config["visualization"]["show_uncertainty"]
                uncertainty_style = config["visualization"]["uncertainty_style"]
                
                if uncertainty_style == "ribbon"
                    # Ribbon style (shaded area)
                    plot!(
                        p,
                        temperatures,
                        values,
                        ribbon=uncertainties,
                        fillalpha=0.2,
                        label=nothing
                    )
                elseif uncertainty_style == "band"
                    # Uncertainty bands (upper and lower curves)
                    plot!(
                        p,
                        temperatures,
                        values + uncertainties,
                        linestyle=:dash,
                        linewidth=1,
                        alpha=0.3,
                        label=nothing
                    )
                    
                    plot!(
                        p,
                        temperatures,
                        values - uncertainties,
                        linestyle=:dash,
                        linewidth=1,
                        alpha=0.3,
                        label=nothing
                    )
                end
            end
        catch e
            @warn "Failed to calculate properties for $species_name: $e"
        end
    end
    
    # Update y-axis label with units if available
    if length(species_names) > 0
        try
            result = calculate_properties(conn, species_names[1], temperature_range, config, step=step)
            units = result["properties"][property_name]["units"]
            ylabel!(p, "$property_name ($(units))")
        catch
            # Keep default label if there's an error
        end
    end
    
    return p
end

"""
    plot_data_source_comparison(conn::DuckDB.DB, species_name::String, 
                              property_name::String, temperature_range::Vector{Float64}, 
                              config::Dict; step::Float64=10.0)

Compare data from different sources for a species.
"""
function plot_data_source_comparison(conn::DuckDB.DB, species_name::String, 
                                   property_name::String, temperature_range::Vector{Float64}, 
                                   config::Dict; step::Float64=10.0)
    # Set up plot theme
    theme = config["visualization"]["default_theme"]
    default_theme = theme == "default" ? :default : Symbol(theme)
    theme!(default_theme)
    
    # Create plot
    p = plot(
        xlabel="Temperature (K)", 
        ylabel="$property_name",
        title="$property_name for $species_name by Source",
        grid=true,
        framestyle=:box
    )
    
    # Get available data sources for this species
    species_info = query_species(conn, species_name, config)
    data_sources = species_info["data_sources"]
    
    # Sort by priority (lowest first, so highest priority appears on top in legend)
    sort!(data_sources, by = s -> s["priority"])
    
    # Generate temperature points
    temperatures = temperature_range[1]:step:temperature_range[2]
    
    # Process each data source
    for source in data_sources
        source_name = source["data_source"]
        temp_min = source["temperature_range"][1]
        temp_max = source["temperature_range"][2]
        
        # Filter temperatures to the range available for this source
        source_temps = filter(t -> temp_min <= t <= temp_max, temperatures)
        
        if isempty(source_temps)
            @warn "No temperature points in range for source: $source_name"
            continue
        end
        
        # Calculate properties
        values = Float64[]
        uncertainties = Float64[]
        units = ""
        
        for temp in source_temps
            try
                # Query with specific source
                data = get_best_thermodynamic_data(conn, species_name, temp)
                
                if data === nothing || data["data_source"] != source_name
                    # Skip if no data or wrong source
                    continue
                end
                
                # Calculate properties
                result = query_properties(conn, species_name, temp, config)
                
                if !haskey(result["properties"], property_name)
                    # Property not available
                    continue
                end
                
                push!(values, result["properties"][property_name]["value"])
                push!(uncertainties, result["properties"][property_name]["uncertainty"])
                units = result["properties"][property_name]["units"]
            catch e
                @warn "Error calculating $property_name at $temp K for source $source_name: $e"
            end
        end
        
        if isempty(values)
            @warn "No valid data points for source: $source_name"
            continue
        end
        
        # Add to plot
        plot!(
            p,
            source_temps[1:length(values)],
            values,
            label=source_name,
            lw=2
        )
        
        # Add uncertainty if enabled
        if config["visualization"]["show_uncertainty"]
            uncertainty_style = config["visualization"]["uncertainty_style"]
            
            if uncertainty_style == "ribbon"
                # Ribbon style (shaded area)
                plot!(
                    p,
                    source_temps[1:length(values)],
                    values,
                    ribbon=uncertainties,
                    fillalpha=0.2,
                    label=nothing
                )
            elseif uncertainty_style == "band"
                # Uncertainty bands (upper and lower curves)
                plot!(
                    p,
                    source_temps[1:length(values)],
                    values + uncertainties,
                    linestyle=:dash,
                    linewidth=1,
                    alpha=0.3,
                    label=nothing
                )
                
                plot!(
                    p,
                    source_temps[1:length(values)],
                    values - uncertainties,
                    linestyle=:dash,
                    linewidth=1,
                    alpha=0.3,
                    label=nothing
                )
            end
        end
        
        # Update y-axis label with units if available
        if !isempty(units)
            ylabel!(p, "$property_name ($(units))")
        end
    end
    
    return p
end

"""
    plot_uncertainty_analysis(conn::DuckDB.DB, species_name::String, 
                            property_name::String, temperature_range::Vector{Float64}, 
                            config::Dict; step::Float64=10.0)

Plot uncertainty analysis for a property.
"""
function plot_uncertainty_analysis(conn::DuckDB.DB, species_name::String, 
                                 property_name::String, temperature_range::Vector{Float64}, 
                                 config::Dict; step::Float64=10.0)
    # Set up plot theme
    theme = config["visualization"]["default_theme"]
    default_theme = theme == "default" ? :default : Symbol(theme)
    theme!(default_theme)
    
    # Calculate properties
    result = calculate_properties(conn, species_name, temperature_range, config, step=step)
    
    temperatures = result["temperatures"]
    values = result["properties"][property_name]["values"]
    uncertainties = result["properties"][property_name]["uncertainties"]
    units = result["properties"][property_name]["units"]
    
    # Calculate relative uncertainty
    relative_uncertainty = 100.0 .* (uncertainties ./ abs.(values))
    
    # Create plots
    p1 = plot(
        temperatures, 
        values, 
        xlabel="Temperature (K)", 
        ylabel="$property_name ($(units))",
        title="$property_name for $species_name",
        ribbon=uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    p2 = plot(
        temperatures, 
        relative_uncertainty, 
        xlabel="Temperature (K)", 
        ylabel="Relative Uncertainty (%)",
        title="Uncertainty Analysis",
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box,
        color=:red
    )
    
    # Combine plots
    p = plot(p1, p2, layout=(2, 1), size=(800, 600))
    
    return p
end