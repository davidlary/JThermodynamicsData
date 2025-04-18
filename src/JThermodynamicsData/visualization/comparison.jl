"""
Functions for visual comparison of thermodynamic data from different sources.
"""

"""
    compare_data_sources_table(conn::DuckDB.DB, species_name::String, 
                             temperature::Float64, config::Dict)

Create a comparison table of thermodynamic data from different sources.
"""
function compare_data_sources_table(conn::DuckDB.DB, species_name::String, 
                                  temperature::Float64, config::Dict)
    # Get available data sources for this species
    species_info = query_species(conn, species_name, config)
    data_sources = species_info["data_sources"]
    
    # Sort by priority (highest first)
    sort!(data_sources, by = s -> -s["priority"])
    
    # Initialize result columns
    sources = String[]
    cp_values = Vector{Union{Float64, Missing}}()
    cp_uncertainties = Vector{Union{Float64, Missing}}()
    h_values = Vector{Union{Float64, Missing}}()
    h_uncertainties = Vector{Union{Float64, Missing}}()
    s_values = Vector{Union{Float64, Missing}}()
    s_uncertainties = Vector{Union{Float64, Missing}}()
    g_values = Vector{Union{Float64, Missing}}()
    g_uncertainties = Vector{Union{Float64, Missing}}()
    reliability_scores = Vector{Union{Float64, Missing}}()
    
    # Process each data source
    for source in data_sources
        source_name = source["data_source"]
        temp_min = source["temperature_range"][1]
        temp_max = source["temperature_range"][2]
        reliability = source["reliability_score"]
        
        # Skip if temperature is outside range
        if temperature < temp_min || temperature > temp_max
            continue
        end
        
        push!(sources, source_name)
        push!(reliability_scores, reliability)
        
        # Query data source
        try
            # Get data directly from the database for this source
            data_query = """
            SELECT 
                td.data_source,
                td.polynomial_type,
                td.data_json,
                td.uncertainty_json
            FROM 
                thermodynamic_data td
            JOIN 
                species s ON td.species_id = s.id
            WHERE 
                s.name = ? AND
                td.data_source = ? AND
                td.temperature_min <= ? AND
                td.temperature_max >= ?
            LIMIT 1
            """
            
            data_result = DuckDB.execute(conn, data_query, 
                                        [species_name, source_name, temperature, temperature])
            data_df = DataFrame(data_result)
            
            if size(data_df, 1) == 0
                # No data found for this source at this temperature
                push!(cp_values, missing)
                push!(cp_uncertainties, missing)
                push!(h_values, missing)
                push!(h_uncertainties, missing)
                push!(s_values, missing)
                push!(s_uncertainties, missing)
                push!(g_values, missing)
                push!(g_uncertainties, missing)
                continue
            end
            
            # Parse data
            row = data_df[1, :]
            polynomial_type = row.polynomial_type
            data_json = JSON.parse(row.data_json)
            
            # Parse uncertainty
            uncertainty_json = Dict()
            if !isempty(row.uncertainty_json)
                try
                    uncertainty_json = JSON.parse(row.uncertainty_json)
                catch e
                    @warn "Failed to parse uncertainty JSON: $e"
                end
            end
            
            # Calculate properties
            if polynomial_type == "nasa7" || polynomial_type == "nasa9"
                # Get coefficients and temperature ranges
                coeffs = data_json["coefficients"]
                temp_ranges = data_json["temperature_ranges"]
                
                # Find applicable range
                range_idx = 0
                for (i, range) in enumerate(temp_ranges)
                    if range[1] <= temperature && temperature <= range[2]
                        range_idx = i
                        break
                    end
                end
                
                if range_idx == 0
                    # This shouldn't happen, but just in case
                    @warn "No applicable temperature range found for $source_name"
                    push!(cp_values, missing)
                    push!(cp_uncertainties, missing)
                    push!(h_values, missing)
                    push!(h_uncertainties, missing)
                    push!(s_values, missing)
                    push!(s_uncertainties, missing)
                    push!(g_values, missing)
                    push!(g_uncertainties, missing)
                    continue
                end
                
                # Get coefficients for this range
                range_coeffs = coeffs[range_idx]
                
                # Get uncertainties if available
                coeffs_uncertainty = zeros(length(range_coeffs))
                if haskey(uncertainty_json, "coefficients_uncertainty")
                    coeff_uncs = uncertainty_json["coefficients_uncertainty"]
                    if length(coeff_uncs) >= range_idx && length(coeff_uncs[range_idx]) == length(range_coeffs)
                        coeffs_uncertainty = coeff_uncs[range_idx]
                    end
                end
                
                # Calculate properties
                if polynomial_type == "nasa7"
                    # Calculate Cp/R
                    cp = calculate_nasa7_cp(range_coeffs, temperature)
                    
                    # Calculate H/RT
                    h = calculate_nasa7_enthalpy(range_coeffs, temperature)
                    
                    # Calculate S/R
                    s = calculate_nasa7_entropy(range_coeffs, temperature)
                else  # nasa9
                    # Calculate Cp/R
                    cp = calculate_nasa9_cp(range_coeffs, temperature)
                    
                    # Calculate H/RT
                    h = calculate_nasa9_enthalpy(range_coeffs, temperature)
                    
                    # Calculate S/R
                    s = calculate_nasa9_entropy(range_coeffs, temperature)
                end
                
                # Calculate uncertainties
                cp_unc = 0.0
                h_unc = 0.0
                s_unc = 0.0
                
                if any(coeffs_uncertainty .> 0)
                    # Linear error propagation (simplified)
                    if polynomial_type == "nasa7"
                        # Cp/R uncertainty
                        cp_unc = sqrt(
                            coeffs_uncertainty[1]^2 +
                            (temperature * coeffs_uncertainty[2])^2 +
                            (temperature^2 * coeffs_uncertainty[3])^2 +
                            (temperature^3 * coeffs_uncertainty[4])^2 +
                            (temperature^4 * coeffs_uncertainty[5])^2
                        )
                        
                        # H/RT uncertainty
                        h_unc = sqrt(
                            coeffs_uncertainty[1]^2 +
                            (temperature/2 * coeffs_uncertainty[2])^2 +
                            (temperature^2/3 * coeffs_uncertainty[3])^2 +
                            (temperature^3/4 * coeffs_uncertainty[4])^2 +
                            (temperature^4/5 * coeffs_uncertainty[5])^2 +
                            (1/temperature * coeffs_uncertainty[6])^2
                        )
                        
                        # S/R uncertainty
                        s_unc = sqrt(
                            (log(temperature) * coeffs_uncertainty[1])^2 +
                            (temperature * coeffs_uncertainty[2])^2 +
                            (temperature^2/2 * coeffs_uncertainty[3])^2 +
                            (temperature^3/3 * coeffs_uncertainty[4])^2 +
                            (temperature^4/4 * coeffs_uncertainty[5])^2 +
                            coeffs_uncertainty[7]^2
                        )
                    else  # nasa9
                        # Simplified uncertainty calculation
                        cp_unc = 0.05 * cp  # 5% uncertainty
                        h_unc = 0.05 * h    # 5% uncertainty
                        s_unc = 0.05 * s    # 5% uncertainty
                    end
                else
                    # Default uncertainties based on reliability
                    cp_unc = (5.0 - reliability) * 0.02 * cp  # 2-10% based on reliability
                    h_unc = (5.0 - reliability) * 0.02 * h    # 2-10% based on reliability
                    s_unc = (5.0 - reliability) * 0.02 * s    # 2-10% based on reliability
                end
                
                # Convert to standard units
                cp_val = cp * R
                cp_unc_val = cp_unc * R
                
                h_val = h * R * temperature / 1000  # kJ/mol
                h_unc_val = h_unc * R * temperature / 1000
                
                s_val = s * R
                s_unc_val = s_unc * R
                
                # Calculate G = H - T*S
                g_val = h_val - temperature * s_val / 1000
                g_unc_val = sqrt(h_unc_val^2 + (temperature * s_unc_val / 1000)^2)
                
                # Store values
                push!(cp_values, cp_val)
                push!(cp_uncertainties, cp_unc_val)
                push!(h_values, h_val)
                push!(h_uncertainties, h_unc_val)
                push!(s_values, s_val)
                push!(s_uncertainties, s_unc_val)
                push!(g_values, g_val)
                push!(g_uncertainties, g_unc_val)
            elseif polynomial_type == "tabular"
                # Get tabular data
                tabular_data = data_json["tabular_data"]
                temps = tabular_data["T"]
                
                # Find exact match or interpolate
                if temperature in temps
                    # Exact match
                    idx = findfirst(t -> t == temperature, temps)
                    
                    # Get properties
                    cp_val = haskey(tabular_data, "Cp") ? tabular_data["Cp"][idx] : missing
                    h_val = haskey(tabular_data, "H") ? tabular_data["H"][idx] : missing
                    s_val = haskey(tabular_data, "S") ? tabular_data["S"][idx] : missing
                    g_val = haskey(tabular_data, "G") ? tabular_data["G"][idx] : missing
                    
                    # Get uncertainties
                    cp_unc_val = haskey(uncertainty_json, "Cp") ? uncertainty_json["Cp"][idx] : 
                                 ismissing(cp_val) ? missing : 0.05 * cp_val
                    
                    h_unc_val = haskey(uncertainty_json, "H") ? uncertainty_json["H"][idx] : 
                               ismissing(h_val) ? missing : 0.05 * h_val
                    
                    s_unc_val = haskey(uncertainty_json, "S") ? uncertainty_json["S"][idx] : 
                               ismissing(s_val) ? missing : 0.05 * s_val
                    
                    g_unc_val = haskey(uncertainty_json, "G") ? uncertainty_json["G"][idx] : 
                               ismissing(g_val) ? missing : 0.05 * g_val
                else
                    # Interpolate
                    # Find closest points
                    lower_idx = findlast(t -> t < temperature, temps)
                    upper_idx = findfirst(t -> t > temperature, temps)
                    
                    if lower_idx === nothing || upper_idx === nothing
                        # Can't interpolate
                        @warn "Can't interpolate for $source_name"
                        push!(cp_values, missing)
                        push!(cp_uncertainties, missing)
                        push!(h_values, missing)
                        push!(h_uncertainties, missing)
                        push!(s_values, missing)
                        push!(s_uncertainties, missing)
                        push!(g_values, missing)
                        push!(g_uncertainties, missing)
                        continue
                    end
                    
                    # Linear interpolation factor
                    t_lower = temps[lower_idx]
                    t_upper = temps[upper_idx]
                    factor = (temperature - t_lower) / (t_upper - t_lower)
                    
                    # Interpolate properties
                    cp_val = haskey(tabular_data, "Cp") ? 
                           tabular_data["Cp"][lower_idx] + factor * (tabular_data["Cp"][upper_idx] - tabular_data["Cp"][lower_idx]) : 
                           missing
                    
                    h_val = haskey(tabular_data, "H") ? 
                           tabular_data["H"][lower_idx] + factor * (tabular_data["H"][upper_idx] - tabular_data["H"][lower_idx]) : 
                           missing
                    
                    s_val = haskey(tabular_data, "S") ? 
                           tabular_data["S"][lower_idx] + factor * (tabular_data["S"][upper_idx] - tabular_data["S"][lower_idx]) : 
                           missing
                    
                    g_val = haskey(tabular_data, "G") ? 
                           tabular_data["G"][lower_idx] + factor * (tabular_data["G"][upper_idx] - tabular_data["G"][lower_idx]) : 
                           missing
                    
                    # Interpolate uncertainties (or use default)
                    cp_unc_val = haskey(uncertainty_json, "Cp") ? 
                               uncertainty_json["Cp"][lower_idx] + factor * (uncertainty_json["Cp"][upper_idx] - uncertainty_json["Cp"][lower_idx]) : 
                               ismissing(cp_val) ? missing : 0.05 * cp_val
                    
                    h_unc_val = haskey(uncertainty_json, "H") ? 
                               uncertainty_json["H"][lower_idx] + factor * (uncertainty_json["H"][upper_idx] - uncertainty_json["H"][lower_idx]) : 
                               ismissing(h_val) ? missing : 0.05 * h_val
                    
                    s_unc_val = haskey(uncertainty_json, "S") ? 
                               uncertainty_json["S"][lower_idx] + factor * (uncertainty_json["S"][upper_idx] - uncertainty_json["S"][lower_idx]) : 
                               ismissing(s_val) ? missing : 0.05 * s_val
                    
                    g_unc_val = haskey(uncertainty_json, "G") ? 
                               uncertainty_json["G"][lower_idx] + factor * (uncertainty_json["G"][upper_idx] - uncertainty_json["G"][lower_idx]) : 
                               ismissing(g_val) ? missing : 0.05 * g_val
                    
                    # Add interpolation uncertainty
                    if !ismissing(cp_unc_val)
                        cp_unc_val += 0.02 * abs(cp_val)  # Add 2% for interpolation
                    end
                    
                    if !ismissing(h_unc_val)
                        h_unc_val += 0.02 * abs(h_val)  # Add 2% for interpolation
                    end
                    
                    if !ismissing(s_unc_val)
                        s_unc_val += 0.02 * abs(s_val)  # Add 2% for interpolation
                    end
                    
                    if !ismissing(g_unc_val)
                        g_unc_val += 0.02 * abs(g_val)  # Add 2% for interpolation
                    end
                end
                
                # Calculate G if not available
                if ismissing(g_val) && !ismissing(h_val) && !ismissing(s_val)
                    g_val = h_val - temperature * s_val / 1000
                    
                    # Propagate uncertainty
                    if !ismissing(h_unc_val) && !ismissing(s_unc_val)
                        g_unc_val = sqrt(h_unc_val^2 + (temperature * s_unc_val / 1000)^2)
                    else
                        g_unc_val = 0.05 * abs(g_val)  # Default 5% uncertainty
                    end
                end
                
                # Store values
                push!(cp_values, cp_val)
                push!(cp_uncertainties, cp_unc_val)
                push!(h_values, h_val)
                push!(h_uncertainties, h_unc_val)
                push!(s_values, s_val)
                push!(s_uncertainties, s_unc_val)
                push!(g_values, g_val)
                push!(g_uncertainties, g_unc_val)
            else
                # Unsupported or unknown polynomial type
                push!(cp_values, missing)
                push!(cp_uncertainties, missing)
                push!(h_values, missing)
                push!(h_uncertainties, missing)
                push!(s_values, missing)
                push!(s_uncertainties, missing)
                push!(g_values, missing)
                push!(g_uncertainties, missing)
            end
        catch e
            @warn "Error processing source $source_name: $e"
            push!(cp_values, missing)
            push!(cp_uncertainties, missing)
            push!(h_values, missing)
            push!(h_uncertainties, missing)
            push!(s_values, missing)
            push!(s_uncertainties, missing)
            push!(g_values, missing)
            push!(g_uncertainties, missing)
        end
    end
    
    # Construct comparison table
    comparison_table = DataFrame(
        Source = sources,
        Reliability = reliability_scores,
        Cp = cp_values,
        Cp_Uncertainty = cp_uncertainties,
        H = h_values,
        H_Uncertainty = h_uncertainties,
        S = s_values,
        S_Uncertainty = s_uncertainties,
        G = g_values,
        G_Uncertainty = g_uncertainties
    )
    
    return comparison_table
end

"""
    compare_to_reference(conn::DuckDB.DB, species_name::String, reference_source::String, 
                       property_name::String, temperature_range::Vector{Float64}, 
                       config::Dict; step::Float64=10.0)

Compare a property from all sources to a reference source.
"""
function compare_to_reference(conn::DuckDB.DB, species_name::String, reference_source::String, 
                            property_name::String, temperature_range::Vector{Float64}, 
                            config::Dict; step::Float64=10.0)
    # Set up plot theme
    theme = config["visualization"]["default_theme"]
    default_theme = theme == "default" ? :default : Symbol(theme)
    theme!(default_theme)
    
    # Get available data sources for this species
    species_info = query_species(conn, species_name, config)
    data_sources = species_info["data_sources"]
    
    # Check if reference source is available
    ref_source_data = nothing
    for source in data_sources
        if source["data_source"] == reference_source
            ref_source_data = source
            break
        end
    end
    
    if ref_source_data === nothing
        error("Reference source $reference_source not available for $species_name")
    end
    
    # Generate temperature points
    temperatures = temperature_range[1]:step:temperature_range[2]
    
    # Filter temperatures to the range available for the reference source
    ref_temp_min = ref_source_data["temperature_range"][1]
    ref_temp_max = ref_source_data["temperature_range"][2]
    ref_temps = filter(t -> ref_temp_min <= t <= ref_temp_max, temperatures)
    
    if isempty(ref_temps)
        error("No temperature points in range for reference source: $reference_source")
    end
    
    # Calculate reference values
    ref_values = Float64[]
    ref_uncertainties = Float64[]
    units = ""
    
    for temp in ref_temps
        try
            # Query with specific source
            data_query = """
            SELECT 
                td.data_source,
                td.polynomial_type,
                td.data_json,
                td.uncertainty_json
            FROM 
                thermodynamic_data td
            JOIN 
                species s ON td.species_id = s.id
            WHERE 
                s.name = ? AND
                td.data_source = ? AND
                td.temperature_min <= ? AND
                td.temperature_max >= ?
            LIMIT 1
            """
            
            data_result = DuckDB.execute(conn, data_query, 
                                       [species_name, reference_source, temp, temp])
            data_df = DataFrame(data_result)
            
            if size(data_df, 1) == 0
                # No data found for this source at this temperature
                continue
            end
            
            # Calculate properties
            result = query_properties(conn, species_name, temp, config)
            
            if !haskey(result["properties"], property_name)
                # Property not available
                continue
            end
            
            push!(ref_values, result["properties"][property_name]["value"])
            push!(ref_uncertainties, result["properties"][property_name]["uncertainty"])
            units = result["properties"][property_name]["units"]
        catch e
            @warn "Error calculating $property_name at $temp K for reference source: $e"
        end
    end
    
    if isempty(ref_values)
        error("Failed to calculate reference values")
    end
    
    # Create plot
    p = plot(
        xlabel="Temperature (K)", 
        ylabel="$property_name Difference ($(units))",
        title="$property_name Deviation from $reference_source",
        grid=true,
        framestyle=:box
    )
    
    # Add reference line at zero
    plot!(
        p,
        ref_temps[1:length(ref_values)],
        zeros(length(ref_values)),
        label="$reference_source",
        lw=2,
        linestyle=:solid,
        color=:black
    )
    
    # Add reference uncertainty band
    plot!(
        p,
        ref_temps[1:length(ref_values)],
        ref_uncertainties,
        fillrange=-ref_uncertainties,
        label="Reference Uncertainty",
        linealpha=0,
        fillalpha=0.2,
        color=:black
    )
    
    # Process each data source (excluding reference)
    for source in data_sources
        source_name = source["data_source"]
        
        # Skip reference source
        if source_name == reference_source
            continue
        end
        
        temp_min = source["temperature_range"][1]
        temp_max = source["temperature_range"][2]
        
        # Filter temperatures to the range available for this source
        source_temps = filter(t -> temp_min <= t <= temp_max, ref_temps[1:length(ref_values)])
        
        if isempty(source_temps)
            @warn "No overlapping temperature points for source: $source_name"
            continue
        end
        
        # Calculate properties and differences
        values = Float64[]
        uncertainties = Float64[]
        
        for temp in source_temps
            try
                # Find reference value for this temperature
                ref_idx = findfirst(t -> t == temp, ref_temps[1:length(ref_values)])
                
                if ref_idx === nothing
                    # No reference value for this temperature
                    continue
                end
                
                ref_value = ref_values[ref_idx]
                
                # Query with specific source
                data_query = """
                SELECT 
                    td.data_source,
                    td.polynomial_type,
                    td.data_json,
                    td.uncertainty_json
                FROM 
                    thermodynamic_data td
                JOIN 
                    species s ON td.species_id = s.id
                WHERE 
                    s.name = ? AND
                    td.data_source = ? AND
                    td.temperature_min <= ? AND
                    td.temperature_max >= ?
                LIMIT 1
                """
                
                data_result = DuckDB.execute(conn, data_query, 
                                           [species_name, source_name, temp, temp])
                data_df = DataFrame(data_result)
                
                if size(data_df, 1) == 0
                    # No data found for this source at this temperature
                    continue
                end
                
                # Calculate properties
                result = query_properties(conn, species_name, temp, config)
                
                if !haskey(result["properties"], property_name)
                    # Property not available
                    continue
                end
                
                source_value = result["properties"][property_name]["value"]
                source_uncertainty = result["properties"][property_name]["uncertainty"]
                
                # Calculate difference from reference
                diff = source_value - ref_value
                
                # Combined uncertainty (root sum of squares)
                combined_uncertainty = sqrt(source_uncertainty^2 + ref_uncertainties[ref_idx]^2)
                
                push!(values, diff)
                push!(uncertainties, combined_uncertainty)
            catch e
                @warn "Error calculating difference at $temp K for source $source_name: $e"
            end
        end
        
        if isempty(values)
            @warn "No valid comparison points for source: $source_name"
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
    end
    
    return p
end