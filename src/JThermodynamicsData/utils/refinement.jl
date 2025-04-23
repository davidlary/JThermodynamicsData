"""
Functions for progressive refinement of thermodynamic data and uncertainty.
"""

"""
    combine_weighted_values(val1::Float64, unc1::Float64, val2::Float64, unc2::Float64, 
                         weight1::Float64, weight2::Float64)

Combine two values with uncertainties using weighted averaging.
Returns the combined value and uncertainty.
"""
function combine_weighted_values(val1::Float64, unc1::Float64, val2::Float64, unc2::Float64, 
                              weight1::Float64, weight2::Float64)
    # Normalize weights
    total_weight = weight1 + weight2
    norm_weight1 = weight1 / total_weight
    norm_weight2 = weight2 / total_weight
    
    # Weighted average for value
    combined_val = norm_weight1 * val1 + norm_weight2 * val2
    
    # Combined uncertainty using error propagation
    # The uncertainty also takes into account the difference between the values
    value_diff_term = (norm_weight1 * norm_weight2 * (val1 - val2)^2)
    uncertainty_term = (norm_weight1 * unc1)^2 + (norm_weight2 * unc2)^2
    
    combined_unc = sqrt(uncertainty_term + value_diff_term)
    
    return combined_val, combined_unc
end

"""
    refine_thermodynamic_data(previous_result::Dict, new_result::Dict, weight_factor::Float64)

Refine thermodynamic data by combining previous result with new data.
Weight factor should be higher for more reliable sources.
Returns the refined data.
"""
function refine_thermodynamic_data(previous_result::Dict, new_result::Dict, weight_factor::Float64)
    # Create a copy of the previous result
    refined_result = deepcopy(previous_result)
    
    # Update source information
    refined_result["data_source"] = new_result["data_source"]
    refined_result["polynomial_type"] = new_result["polynomial_type"]
    
    # Check if this is an atomic or ionic species
    formula = get(new_result, "formula", "")
    is_atomic_or_ion = (!isempty(formula) && 
                        (any(f -> f == formula, ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
                                               "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
                                               "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", 
                                               "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
                                               "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", 
                                               "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"]) || 
                         endswith(formula, "+") || endswith(formula, "-")))
    
    # For atomic species and ions with experimental data, use highest priority data directly
    # if weight_factor is close to 1.0, it means this is a high-priority source
    if is_atomic_or_ion && !startswith(new_result["data_source"], "THEORETICAL") && weight_factor > 0.9
        @info "Using highest priority source $(new_result["data_source"]) directly for atomic/ionic species $formula"
        # Take the new result values directly for all properties
        for prop in ["Cp", "H", "S", "G"]
            if haskey(new_result["properties"], prop) && haskey(refined_result["properties"], prop)
                refined_result["properties"][prop]["value"] = new_result["properties"][prop]["value"]
                refined_result["properties"][prop]["uncertainty"] = new_result["properties"][prop]["uncertainty"]
            end
        end
        
        # Remove theoretical method info if using experimental data
        if haskey(refined_result, "theoretical_method")
            delete!(refined_result, "theoretical_method")
        end
        
        return refined_result
    end
    
    # Remove theoretical method info if using experimental data
    if haskey(refined_result, "theoretical_method")
        delete!(refined_result, "theoretical_method")
    end
    
    # Calculate weights - higher weight for new data from more reliable source
    previous_weight = 1.0 - weight_factor  # Lower weight for previous data
    new_weight = weight_factor  # Higher weight for new data
    
    # Update each primary property
    for prop in ["Cp", "H", "S"]
        if haskey(new_result["properties"], prop) && haskey(previous_result["properties"], prop)
            # Get current and new values
            current_val = previous_result["properties"][prop]["value"]
            current_unc = previous_result["properties"][prop]["uncertainty"]
            
            new_val = new_result["properties"][prop]["value"]
            new_unc = new_result["properties"][prop]["uncertainty"]
            
            # Combine values and uncertainties using weighted averaging
            combined_val, combined_unc = combine_weighted_values(
                current_val, current_unc, new_val, new_unc, previous_weight, new_weight
            )
            
            # Update refined result
            refined_result["properties"][prop]["value"] = combined_val
            refined_result["properties"][prop]["uncertainty"] = combined_unc
        end
    end
    
    # Get temperature for recalculating G
    temperature = get(new_result, "temperature", get(previous_result, "temperature", 298.15))
    
    # Calculate G from H and S consistently to maintain thermodynamic relationship
    if haskey(refined_result["properties"], "H") && haskey(refined_result["properties"], "S")
        h_val = refined_result["properties"]["H"]["value"]
        h_unc = refined_result["properties"]["H"]["uncertainty"]
        s_val = refined_result["properties"]["S"]["value"]
        s_unc = refined_result["properties"]["S"]["uncertainty"]
        
        # Calculate G = H - T*S consistently
        g_val = h_val - temperature * s_val / 1000
        g_unc = sqrt(h_unc^2 + (temperature * s_unc / 1000)^2)
        
        # Update G in the refined result
        refined_result["properties"]["G"]["value"] = g_val
        refined_result["properties"]["G"]["uncertainty"] = g_unc
    else
        # If H or S is missing, update G directly
        if haskey(new_result["properties"], "G") && haskey(previous_result["properties"], "G")
            current_val = previous_result["properties"]["G"]["value"]
            current_unc = previous_result["properties"]["G"]["uncertainty"]
            
            new_val = new_result["properties"]["G"]["value"]
            new_unc = new_result["properties"]["G"]["uncertainty"]
            
            # Combine values and uncertainties using weighted averaging
            combined_val, combined_unc = combine_weighted_values(
                current_val, current_unc, new_val, new_unc, previous_weight, new_weight
            )
            
            # Update refined result
            refined_result["properties"]["G"]["value"] = combined_val
            refined_result["properties"]["G"]["uncertainty"] = combined_unc
        end
    end
    
    return refined_result
end

"""
    get_theoretical_estimate(formula::String, temperature::Float64, config::Dict)

Generate theoretical estimates of thermodynamic properties.
"""
function get_theoretical_estimate(formula::String, temperature::Float64, metadata::Dict, config::Dict)
    # Try different theoretical methods in order of priority
    methods = config["theoretical_calculation"]["methods"]
    
    # Sort by priority
    sort!(methods, by = m -> m["priority"])
    
    # Try each method until one works
    theoretical_result = nothing
    method_used = ""
    
    for method in methods
        if !method["enabled"]
            continue
        end
        
        method_name = method["name"]
        
        try
            if method_name == "group_contribution"
                theoretical_result = estimate_properties_group_contribution(formula, temperature)
                method_used = "group_contribution"
                break
            elseif method_name == "statistical_thermodynamics"
                # Try to get molecular data from metadata
                molecular_data = Dict{String, Any}()
                
                if haskey(metadata, "molecular_data")
                    molecular_data = metadata["molecular_data"]
                else
                    # Estimate molecular properties
                    molecular_data = estimate_molecular_properties(formula)
                end
                
                theoretical_result = estimate_properties_statistical_thermodynamics(
                    formula, molecular_data, temperature
                )
                method_used = "statistical_thermodynamics"
                break
            elseif method_name == "machine_learning"
                theoretical_result = estimate_properties_machine_learning(formula, temperature, config)
                method_used = "machine_learning"
                break
            elseif method_name == "quantum_chemistry"
                theoretical_result = estimate_properties_quantum_chemistry(formula, temperature, config)
                method_used = "quantum_chemistry"
                break
            end
        catch e
            @warn "Failed to calculate properties using $method_name: $e"
        end
    end
    
    if theoretical_result === nothing
        @warn "All theoretical methods failed for formula $formula at $temperature K"
        
        # Create empty theoretical result with high uncertainty
        theoretical_result = Dict(
            "temperature" => temperature,
            "Cp" => 0.0,
            "Cp_uncertainty" => 1000.0,
            "H" => 0.0,
            "H_uncertainty" => 1000.0,
            "S" => 0.0,
            "S_uncertainty" => 1000.0,
            "G" => 0.0,
            "G_uncertainty" => 1000.0,
            "method" => "none"
        )
        method_used = "none"
    end
    
    # Create result structure
    result = Dict(
        "data_source" => "THEORETICAL_$(uppercase(method_used))",
        "polynomial_type" => "theoretical",
        "properties" => Dict(
            "Cp" => Dict(
                "value" => theoretical_result["Cp"],
                "uncertainty" => theoretical_result["Cp_uncertainty"],
                "units" => "J/mol/K"
            ),
            "H" => Dict(
                "value" => theoretical_result["H"],
                "uncertainty" => theoretical_result["H_uncertainty"],
                "units" => "kJ/mol"
            ),
            "S" => Dict(
                "value" => theoretical_result["S"],
                "uncertainty" => theoretical_result["S_uncertainty"],
                "units" => "J/mol/K"
            ),
            "G" => Dict(
                "value" => theoretical_result["G"],
                "uncertainty" => theoretical_result["G_uncertainty"],
                "units" => "kJ/mol"
            )
        ),
        "theoretical_method" => method_used
    )
    
    return result
end

"""
    progressively_refine_thermodynamic_data(conn::DuckDB.DB, species_name::String, 
                                        temperature::Float64, config::Dict)

Progressively refine thermodynamic data by traversing the hierarchical data sources.
This version ensures thorough traversal of all data sources for maximum accuracy.
"""
function progressively_refine_thermodynamic_data(conn::DuckDB.DB, species_name::String, 
                                             temperature::Float64, config::Dict)
    # Get species information
    species_query = """
    SELECT 
        s.id,
        s.formula,
        s.molecular_weight,
        s.metadata_json
    FROM 
        species s
    WHERE 
        s.name = ?
    """
    
    species_result = DuckDB.execute(conn, species_query, [species_name])
    species_df = DataFrame(species_result)
    
    if size(species_df, 1) == 0
        error("Species not found: $species_name")
    end
    
    species_id = species_df[1, :id]
    formula = species_df[1, :formula]
    mw = species_df[1, :molecular_weight]
    
    # Parse metadata
    metadata = Dict{String, Any}()
    if !isempty(species_df[1, :metadata_json]) && species_df[1, :metadata_json] !== missing
        try
            metadata = JSON.parse(species_df[1, :metadata_json])
        catch e
            @warn "Failed to parse metadata JSON: $e"
        end
    end
    
    # STEP 1: Start with multiple theoretical calculations
    @info "Starting with theoretical calculations for $species_name at $temperature K"
    
    # Try all enabled theoretical methods
    all_theoretical_results = Dict[]
    methods = config["theoretical_calculation"]["methods"]
    
    # Try each method and collect all successful results
    for method in methods
        if !method["enabled"]
            continue
        end
        
        method_name = method["name"]
        
        try
            if method_name == "group_contribution"
                theo_result = estimate_properties_group_contribution(formula, temperature)
                
                # Create method result
                method_result = Dict(
                    "species_name" => species_name,
                    "formula" => formula,
                    "temperature" => temperature,
                    "data_source" => "THEORETICAL_GROUP_CONTRIBUTION",
                    "polynomial_type" => "theoretical",
                    "properties" => Dict(
                        "Cp" => Dict(
                            "value" => theo_result["Cp"],
                            "uncertainty" => theo_result["Cp_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "H" => Dict(
                            "value" => theo_result["H"],
                            "uncertainty" => theo_result["H_uncertainty"],
                            "units" => "kJ/mol"
                        ),
                        "S" => Dict(
                            "value" => theo_result["S"],
                            "uncertainty" => theo_result["S_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "G" => Dict(
                            "value" => theo_result["G"],
                            "uncertainty" => theo_result["G_uncertainty"],
                            "units" => "kJ/mol"
                        )
                    ),
                    "theoretical_method" => method_name
                )
                
                push!(all_theoretical_results, method_result)
            elseif method_name == "statistical_thermodynamics"
                # Try to get molecular data from metadata
                molecular_data = Dict{String, Any}()
                
                if haskey(metadata, "molecular_data")
                    molecular_data = metadata["molecular_data"]
                else
                    # Estimate molecular properties
                    molecular_data = estimate_molecular_properties(formula)
                end
                
                theo_result = estimate_properties_statistical_thermodynamics(
                    formula, molecular_data, temperature
                )
                
                # Create method result
                method_result = Dict(
                    "species_name" => species_name,
                    "formula" => formula,
                    "temperature" => temperature,
                    "data_source" => "THEORETICAL_STATISTICAL_THERMO",
                    "polynomial_type" => "theoretical",
                    "properties" => Dict(
                        "Cp" => Dict(
                            "value" => theo_result["Cp"],
                            "uncertainty" => theo_result["Cp_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "H" => Dict(
                            "value" => theo_result["H"],
                            "uncertainty" => theo_result["H_uncertainty"],
                            "units" => "kJ/mol"
                        ),
                        "S" => Dict(
                            "value" => theo_result["S"],
                            "uncertainty" => theo_result["S_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "G" => Dict(
                            "value" => theo_result["G"],
                            "uncertainty" => theo_result["G_uncertainty"],
                            "units" => "kJ/mol"
                        )
                    ),
                    "theoretical_method" => method_name
                )
                
                push!(all_theoretical_results, method_result)
            elseif method_name == "machine_learning"
                theo_result = estimate_properties_machine_learning(formula, temperature, config)
                
                # Create method result
                method_result = Dict(
                    "species_name" => species_name,
                    "formula" => formula,
                    "temperature" => temperature,
                    "data_source" => "THEORETICAL_MACHINE_LEARNING",
                    "polynomial_type" => "theoretical",
                    "properties" => Dict(
                        "Cp" => Dict(
                            "value" => theo_result["Cp"],
                            "uncertainty" => theo_result["Cp_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "H" => Dict(
                            "value" => theo_result["H"],
                            "uncertainty" => theo_result["H_uncertainty"],
                            "units" => "kJ/mol"
                        ),
                        "S" => Dict(
                            "value" => theo_result["S"],
                            "uncertainty" => theo_result["S_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "G" => Dict(
                            "value" => theo_result["G"],
                            "uncertainty" => theo_result["G_uncertainty"],
                            "units" => "kJ/mol"
                        )
                    ),
                    "theoretical_method" => method_name
                )
                
                push!(all_theoretical_results, method_result)
            elseif method_name == "quantum_chemistry"
                theo_result = estimate_properties_quantum_chemistry(formula, temperature, config)
                
                # Create method result
                method_result = Dict(
                    "species_name" => species_name,
                    "formula" => formula,
                    "temperature" => temperature,
                    "data_source" => "THEORETICAL_QUANTUM_CHEMISTRY",
                    "polynomial_type" => "theoretical",
                    "properties" => Dict(
                        "Cp" => Dict(
                            "value" => theo_result["Cp"],
                            "uncertainty" => theo_result["Cp_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "H" => Dict(
                            "value" => theo_result["H"],
                            "uncertainty" => theo_result["H_uncertainty"],
                            "units" => "kJ/mol"
                        ),
                        "S" => Dict(
                            "value" => theo_result["S"],
                            "uncertainty" => theo_result["S_uncertainty"],
                            "units" => "J/mol/K"
                        ),
                        "G" => Dict(
                            "value" => theo_result["G"],
                            "uncertainty" => theo_result["G_uncertainty"],
                            "units" => "kJ/mol"
                        )
                    ),
                    "theoretical_method" => method_name
                )
                
                push!(all_theoretical_results, method_result)
            end
        catch e
            @warn "Failed to calculate properties using $method_name: $e"
        end
    end
    
    # If we have multiple theoretical results, use them to improve uncertainty estimates
    result = Dict()
    if length(all_theoretical_results) > 0
        # Use the first method as a starting point
        result = deepcopy(all_theoretical_results[1])
        
        # If we have multiple methods, refine using their differences
        if length(all_theoretical_results) > 1
            @info "Using $(length(all_theoretical_results)) theoretical methods to refine uncertainty"
            
            # For each property, calculate mean and standard deviation across methods
            for prop in ["Cp", "H", "S", "G"]
                values = Float64[]
                
                for method_result in all_theoretical_results
                    if haskey(method_result["properties"], prop)
                        push!(values, method_result["properties"][prop]["value"])
                    end
                end
                
                if length(values) > 1
                    # Calculate mean and standard deviation
                    mean_val = mean(values)
                    std_val = std(values)
                    
                    # Update result with mean value
                    result["properties"][prop]["value"] = mean_val
                    
                    # Use standard deviation between methods as uncertainty
                    # But don't reduce uncertainty below the minimum from individual methods
                    min_uncertainty = minimum([r["properties"][prop]["uncertainty"] for r in all_theoretical_results 
                                             if haskey(r["properties"], prop)])
                    
                    # Set uncertainty to max of std between methods or minimum method uncertainty
                    result["properties"][prop]["uncertainty"] = max(std_val, min_uncertainty)
                end
            end
            
            # Set source to indicate multiple methods
            result["data_source"] = "THEORETICAL_ENSEMBLE"
            result["theoretical_method"] = "ensemble"
        end
    else
        # No theoretical methods worked, create fallback with high uncertainty
        @warn "All theoretical methods failed for $species_name at $temperature K"
        
        result = Dict(
            "species_name" => species_name,
            "formula" => formula,
            "temperature" => temperature,
            "data_source" => "THEORETICAL_FALLBACK",
            "polynomial_type" => "theoretical",
            "properties" => Dict(
                "Cp" => Dict(
                    "value" => 0.0,
                    "uncertainty" => 1000.0,
                    "units" => "J/mol/K"
                ),
                "H" => Dict(
                    "value" => 0.0,
                    "uncertainty" => 1000.0,
                    "units" => "kJ/mol"
                ),
                "S" => Dict(
                    "value" => 0.0,
                    "uncertainty" => 1000.0,
                    "units" => "J/mol/K"
                ),
                "G" => Dict(
                    "value" => 0.0,
                    "uncertainty" => 1000.0,
                    "units" => "kJ/mol"
                )
            ),
            "theoretical_method" => "fallback"
        )
    end
    
    # Store all theoretical sources for the record
    all_sources = all_theoretical_results
    
    # STEP 2: Traverse the hierarchy of data sources
    @info "Traversing hierarchy of data sources for $species_name at $temperature K"
    
    # Get ALL data sources for this species, not just ones that perfectly match the temperature
    # We'll check each source to see if it can be used for this temperature
    local all_thermo_data_df = DataFrame()
    try
        # First get sources with temperature ranges that include our target temperature
        all_thermo_data_df = get_all_thermodynamic_data_for_species(conn, species_name)
    catch e
        @warn "Error getting thermodynamic data: $e"
        all_thermo_data_df = DataFrame()
    end
    
    # Create a list of all data sources from the hierarchy
    # This ensures we check for the presence of every source in the hierarchy
    hierarchy_sources = [
        "GRI-MECH",         # Priority 1
        "CHEMKIN",          # Priority 2
        "NASA-CEA",         # Priority 3
        "JANAF",            # Priority 4
        "THERMOML",         # Priority 5
        "TDE",              # Priority 6
        "BURCAT",           # Priority 7
        "ATCT"              # Priority 8
    ]
    
    # Filter the data sources that are applicable for this temperature
    applicable_data_df = filter(row -> 
        row.temperature_min <= temperature && temperature <= row.temperature_max, 
        all_thermo_data_df
    )
    
    # Sort by priority (high to low)
    sort!(applicable_data_df, :priority, rev=true)
    
    # Record sources found at each priority level
    sources_by_priority = Dict{Int, Vector{String}}()
    for row in eachrow(applicable_data_df)
        priority = row.priority
        source = row.data_source
        
        if !haskey(sources_by_priority, priority)
            sources_by_priority[priority] = String[]
        end
        
        push!(sources_by_priority[priority], source)
    end
    
    # Log the sources found at each priority level
    for priority in sort(collect(keys(sources_by_priority)), rev=true)
        sources = sources_by_priority[priority]
        @info "Found $(length(sources)) source(s) at priority level $priority: $(join(sources, ", "))"
    end
    
    # Check which sources from the hierarchy are missing
    sources_found = Set([row.data_source for row in eachrow(applicable_data_df)])
    missing_sources = setdiff(Set(hierarchy_sources), sources_found)
    
    if !isempty(missing_sources)
        @info "Missing sources for $species_name at $temperature K: $(join(collect(missing_sources), ", "))"
    end
    
    # Track if any sources have been used
    experimental_sources_used = false
    
    if size(applicable_data_df, 1) > 0
        # Process each data source to progressively refine the data and uncertainty
        for (i, row) in enumerate(eachrow(applicable_data_df))
            source_name = row.data_source
            priority = row.priority
            @info "Processing data source: $source_name ($(i)/$(size(applicable_data_df, 1)))"
            
            # Parse data from this source
            data_json = JSON.parse(row.data_json)
            
            # Parse uncertainty if available
            uncertainty_json = Dict()
            if !isempty(row.uncertainty_json)
                try
                    uncertainty_json = JSON.parse(row.uncertainty_json)
                catch e
                    @warn "Failed to parse uncertainty JSON: $e"
                end
            end
            
            # Calculate properties from this source
            source_result = Dict(
                "species_name" => species_name,
                "formula" => formula,
                "temperature" => temperature,
                "data_source" => source_name,
                "polynomial_type" => row.polynomial_type,
                "properties" => Dict()
            )
            
            # Calculate properties based on polynomial type
            if row.polynomial_type == "nasa7" || row.polynomial_type == "nasa9"
                # Get coefficients and temperature ranges
                coeffs = data_json["coefficients"]
                temp_ranges = data_json["temperature_ranges"]
                
                # Find applicable range
                range_idx = 0
                for (idx, range) in enumerate(temp_ranges)
                    if range[1] <= temperature && temperature <= range[2]
                        range_idx = idx
                        break
                    end
                end
                
                if range_idx == 0
                    @warn "No applicable temperature range found for $source_name"
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
                if row.polynomial_type == "nasa7"
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
                    if row.polynomial_type == "nasa7"
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
                    reliability = row.reliability_score
                    
                    # Calculate uncertainty based on reliability score and priority
                    # Higher priority sources (closer to 8) get lower uncertainty
                    priority_factor = max(0.5, 1.0 - (priority / 10.0))  # 0.5 to 0.9 for priorities 1-8
                    
                    cp_unc = (5.0 - reliability) * 0.02 * cp * priority_factor  # 1-5% based on reliability and priority
                    h_unc = (5.0 - reliability) * 0.02 * h * priority_factor    # 1-5% based on reliability and priority
                    s_unc = (5.0 - reliability) * 0.02 * s * priority_factor    # 1-5% based on reliability and priority
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
                
                # Store properties
                source_result["properties"]["Cp"] = Dict(
                    "value" => cp_val,
                    "uncertainty" => cp_unc_val,
                    "units" => "J/mol/K"
                )
                
                source_result["properties"]["H"] = Dict(
                    "value" => h_val,
                    "uncertainty" => h_unc_val,
                    "units" => "kJ/mol"
                )
                
                source_result["properties"]["S"] = Dict(
                    "value" => s_val,
                    "uncertainty" => s_unc_val,
                    "units" => "J/mol/K"
                )
                
                source_result["properties"]["G"] = Dict(
                    "value" => g_val,
                    "uncertainty" => g_unc_val,
                    "units" => "kJ/mol"
                )
                
                # Mark that we've used an experimental source
                experimental_sources_used = true
            elseif row.polynomial_type == "tabular"
                # Process tabular data (JANAF, etc.)
                # For JANAF and other tabular data, we need to interpolate between data points
                
                if haskey(data_json, "data_points")
                    data_points = data_json["data_points"]
                    
                    # Sort data points by temperature
                    sort!(data_points, by = p -> p["temperature"])
                    
                    # Find closest points for interpolation
                    lower_idx = 0
                    upper_idx = 0
                    
                    for (idx, point) in enumerate(data_points)
                        point_temp = point["temperature"]
                        
                        if point_temp <= temperature
                            lower_idx = idx
                        end
                        
                        if point_temp >= temperature && upper_idx == 0
                            upper_idx = idx
                        end
                    end
                    
                    # If temperature exactly matches a data point
                    if lower_idx > 0 && data_points[lower_idx]["temperature"] == temperature
                        # Directly use this data point
                        point = data_points[lower_idx]
                        
                        # Extract properties
                        cp_val = get(point, "Cp", 0.0)
                        h_val = get(point, "H", 0.0)
                        s_val = get(point, "S", 0.0)
                        g_val = get(point, "G", h_val - temperature * s_val / 1000)
                        
                        # Get uncertainties if available
                        cp_unc_val = get(point, "Cp_uncertainty", 0.05 * cp_val)
                        h_unc_val = get(point, "H_uncertainty", 0.05 * h_val)
                        s_unc_val = get(point, "S_uncertainty", 0.05 * s_val)
                        g_unc_val = sqrt(h_unc_val^2 + (temperature * s_unc_val / 1000)^2)
                        
                        # Store properties
                        source_result["properties"]["Cp"] = Dict(
                            "value" => cp_val,
                            "uncertainty" => cp_unc_val,
                            "units" => "J/mol/K"
                        )
                        
                        source_result["properties"]["H"] = Dict(
                            "value" => h_val,
                            "uncertainty" => h_unc_val,
                            "units" => "kJ/mol"
                        )
                        
                        source_result["properties"]["S"] = Dict(
                            "value" => s_val,
                            "uncertainty" => s_unc_val,
                            "units" => "J/mol/K"
                        )
                        
                        source_result["properties"]["G"] = Dict(
                            "value" => g_val,
                            "uncertainty" => g_unc_val,
                            "units" => "kJ/mol"
                        )
                        
                        # Mark that we've used an experimental source
                        experimental_sources_used = true
                    elseif lower_idx > 0 && upper_idx > 0 && lower_idx != upper_idx
                        # Interpolate between points
                        lower_point = data_points[lower_idx]
                        upper_point = data_points[upper_idx]
                        
                        lower_temp = lower_point["temperature"]
                        upper_temp = upper_point["temperature"]
                        
                        # Linear interpolation factor
                        factor = (temperature - lower_temp) / (upper_temp - lower_temp)
                        
                        # Interpolate properties
                        cp_val = lower_point["Cp"] + factor * (upper_point["Cp"] - lower_point["Cp"])
                        h_val = lower_point["H"] + factor * (upper_point["H"] - lower_point["H"])
                        s_val = lower_point["S"] + factor * (upper_point["S"] - lower_point["S"])
                        g_val = h_val - temperature * s_val / 1000
                        
                        # Interpolate or estimate uncertainties
                        cp_unc_val = get(lower_point, "Cp_uncertainty", 0.05 * cp_val)
                        h_unc_val = get(lower_point, "H_uncertainty", 0.05 * h_val)
                        s_unc_val = get(lower_point, "S_uncertainty", 0.05 * s_val)
                        g_unc_val = sqrt(h_unc_val^2 + (temperature * s_unc_val / 1000)^2)
                        
                        # Store properties
                        source_result["properties"]["Cp"] = Dict(
                            "value" => cp_val,
                            "uncertainty" => cp_unc_val,
                            "units" => "J/mol/K"
                        )
                        
                        source_result["properties"]["H"] = Dict(
                            "value" => h_val,
                            "uncertainty" => h_unc_val,
                            "units" => "kJ/mol"
                        )
                        
                        source_result["properties"]["S"] = Dict(
                            "value" => s_val,
                            "uncertainty" => s_unc_val,
                            "units" => "J/mol/K"
                        )
                        
                        source_result["properties"]["G"] = Dict(
                            "value" => g_val,
                            "uncertainty" => g_unc_val,
                            "units" => "kJ/mol"
                        )
                        
                        # Mark that we've used an experimental source
                        experimental_sources_used = true
                    else
                        @warn "Cannot interpolate tabular data for $source_name at $temperature K"
                        continue
                    end
                else
                    @warn "Tabular data source $source_name has no data points"
                    continue
                end
            else
                @warn "Unsupported polynomial type: $(row.polynomial_type)"
                continue
            end
            
            # Calculate weight factor for this source
            # Weight increases with priority, so highest priority (ATcT, Priority 8) gets highest weight
            weight_factor = 0.5 + 0.5 * (priority / 8.0)  # Increases from 0.5 to 1.0 based on priority
            
            # Store this source result
            push!(all_sources, source_result)
            
            # Refine the result by combining with this source
            result = refine_thermodynamic_data(result, source_result, weight_factor)
        end
    end
    
    # Check if we used any experimental data sources
    if !experimental_sources_used
        @warn "No experimental data sources used for $species_name at $temperature K, using theoretical estimates only"
    end
    
    # Store the result in cache for future reference
    store_in_cache(conn, species_id, temperature, temperature, "properties", result)
    
    # Store data source usage information for documentation
    if !haskey(result, "source_hierarchy")
        result["source_hierarchy"] = Dict(
            "theoretical_methods" => [r["theoretical_method"] for r in all_theoretical_results],
            "experimental_sources" => [r["data_source"] for r in all_sources if startswith(r["data_source"], "THEORETICAL") == false],
            "missing_sources" => collect(missing_sources)
        )
    end
    
    # Return the result and all sources used
    return result, all_sources
end

"""
    get_all_thermodynamic_data_for_species(conn::DuckDB.DB, species_name::String)

Get all available thermodynamic data for a species across all sources.
This ensures we thoroughly check every possible data source for the species.
"""
function get_all_thermodynamic_data_for_species(conn::DuckDB.DB, species_name::String)
    # Query all thermodynamic data for the species
    query = """
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
        ds.priority,
        td.date_modified
    FROM 
        species s
    JOIN 
        thermodynamic_data td ON s.id = td.species_id
    JOIN 
        data_sources ds ON td.data_source = ds.name
    WHERE 
        s.name = ?
    ORDER BY 
        ds.priority DESC, 
        td.reliability_score DESC,
        td.date_modified DESC
    """
    
    result = DuckDB.execute(conn, query, [species_name])
    return DataFrame(result)
end

"""
    generate_source_usage_report(conn::DuckDB.DB, species_list::Vector{String}, temperature::Float64, config::Dict)

Generate a report on data source usage for a list of species.
This helps identify which sources are used most frequently and which species lack experimental data.
"""
function generate_source_usage_report(conn::DuckDB.DB, species_list::Vector{String}, temperature::Float64, config::Dict)
    # Initialize counters
    source_usage = Dict{String, Int}()
    theoretical_only = String[]
    missing_data = String[]
    
    # Process each species
    for species_name in species_list
        try
            # Query properties for this species
            result, sources = progressively_refine_thermodynamic_data(conn, species_name, temperature, config)
            
            # Count sources used
            experimental_sources = filter(s -> !startswith(s["data_source"], "THEORETICAL"), sources)
            
            if isempty(experimental_sources)
                push!(theoretical_only, species_name)
            end
            
            # Update source usage counters
            for source in sources
                source_name = source["data_source"]
                source_usage[source_name] = get(source_usage, source_name, 0) + 1
            end
        catch e
            @warn "Error processing $species_name: $e"
            push!(missing_data, species_name)
        end
    end
    
    # Create report
    report = Dict(
        "temperature" => temperature,
        "total_species" => length(species_list),
        "processed_species" => length(species_list) - length(missing_data),
        "theoretical_only_count" => length(theoretical_only),
        "theoretical_only_species" => theoretical_only,
        "missing_data" => missing_data,
        "source_usage" => source_usage
    )
    
    return report
end

"""
    create_markdown_documentation(species_name::String, result::Dict, all_sources::Vector{Dict}, output_dir::String)

Create detailed Markdown documentation for a species detailing all data sources used.
"""
function create_markdown_documentation(species_name::String, result::Dict, all_sources::Vector{Dict}, output_dir::String)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Get formula
    formula = result["formula"]
    temperature = result["temperature"]
    
    # Sort sources by priority
    # First theoretical, then experimental by data source priority
    sorted_sources = sort(all_sources, by = s -> startswith(s["data_source"], "THEORETICAL") ? 0 : 
                          (haskey(s, "priority") ? s["priority"] : 9))
    
    # Count theoretical vs experimental sources
    theoretical_count = count(s -> startswith(s["data_source"], "THEORETICAL"), sorted_sources)
    experimental_count = length(sorted_sources) - theoretical_count
    
    # Create markdown file
    markdown_file = joinpath(output_dir, "$(species_name)_sources.md")
    open(markdown_file, "w") do io
        write(io, "# Thermodynamic Data Sources for $species_name ($formula)\n\n")
        
        # Summary
        write(io, "## Summary\n\n")
        write(io, "- Temperature: $(temperature) K\n")
        write(io, "- Theoretical methods used: $theoretical_count\n")
        write(io, "- Experimental sources used: $experimental_count\n\n")
        
        # Final properties
        write(io, "## Final Properties\n\n")
        write(io, "| Property | Value | Uncertainty | Units |\n")
        write(io, "|----------|-------|------------|-------|\n")
        
        for prop in ["Cp", "H", "S", "G"]
            if haskey(result["properties"], prop)
                value = result["properties"][prop]["value"]
                uncertainty = result["properties"][prop]["uncertainty"]
                units = result["properties"][prop]["units"]
                write(io, "| $prop | $value | $uncertainty | $units |\n")
            end
        end
        
        write(io, "\n")
        
        # Source details
        write(io, "## Data Sources\n\n")
        write(io, "Listed in order of application, from theoretical baselines to highest priority experimental:\n\n")
        
        # Create table of sources
        write(io, "| Source | Type | Properties |\n")
        write(io, "|--------|------|------------|\n")
        
        for source in sorted_sources
            source_name = source["data_source"]
            source_type = source["polynomial_type"]
            
            # Summarize properties
            props = []
            for prop in ["Cp", "H", "S", "G"]
                if haskey(source["properties"], prop)
                    value = round(source["properties"][prop]["value"], digits=2)
                    push!(props, "$prop=$value")
                end
            end
            
            properties_str = join(props, ", ")
            write(io, "| $source_name | $source_type | $properties_str |\n")
        end
        
        # Source hierarchy info if available
        if haskey(result, "source_hierarchy")
            write(io, "\n## Hierarchy Information\n\n")
            
            if haskey(result["source_hierarchy"], "missing_sources") && !isempty(result["source_hierarchy"]["missing_sources"])
                write(io, "### Missing Sources\n\n")
                for source in result["source_hierarchy"]["missing_sources"]
                    write(io, "- $source\n")
                end
                write(io, "\n")
            end
        end
        
        # Generation info
        write(io, "\n---\n")
        write(io, "*Generated on $(Dates.now()) using JThermodynamicsData hierarchical calculation system*\n")
    end
    
    return markdown_file
end

"""
    create_json_documentation(species_name::String, result::Dict, all_sources::Vector{Dict}, output_dir::String)

Create detailed JSON documentation for a species detailing all data sources used.
"""
function create_json_documentation(species_name::String, result::Dict, all_sources::Vector{Dict}, output_dir::String)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Create documentation JSON
    doc = Dict(
        "species_name" => species_name,
        "formula" => result["formula"],
        "temperature" => result["temperature"],
        "final_properties" => result["properties"],
        "all_sources" => all_sources,
        "generation_date" => string(Dates.now())
    )
    
    # Add source hierarchy info if available
    if haskey(result, "source_hierarchy")
        doc["source_hierarchy"] = result["source_hierarchy"]
    end
    
    # Write JSON file
    json_file = joinpath(output_dir, "$(species_name)_sources.json")
    open(json_file, "w") do io
        write(io, JSON.json(doc, 2))  # Pretty-printed with 2-space indentation
    end
    
    return json_file
end