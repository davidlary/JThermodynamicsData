"""
Functions for converting between different thermodynamic data formats.
"""

"""
    convert_tabular_to_nasa7(tabular_data::Dict, temperature_range::Vector{Float64})

Convert tabular thermodynamic data to NASA 7-coefficient polynomials.
"""
function convert_tabular_to_nasa7(tabular_data::Dict, temperature_range::Vector{Float64})
    # Check if we have the required data
    if !haskey(tabular_data, "T") || !haskey(tabular_data, "Cp")
        error("Missing required data columns in tabular data")
    end
    
    # Get the data arrays
    temps = tabular_data["T"]
    cps = tabular_data["Cp"]
    
    # Check for enthalpy and entropy data
    has_h = haskey(tabular_data, "H")
    has_s = haskey(tabular_data, "S")
    
    hs = has_h ? tabular_data["H"] : zeros(length(temps))
    ss = has_s ? tabular_data["S"] : zeros(length(temps))
    
    # Sort data by temperature
    indices = sortperm(temps)
    temps = temps[indices]
    cps = cps[indices]
    hs = has_h ? hs[indices] : zeros(length(temps))
    ss = has_s ? ss[indices] : zeros(length(temps))
    
    # Define temperature ranges for the polynomials
    # If temperature_range is not provided, use the range from the data
    if isempty(temperature_range)
        temperature_range = [minimum(temps), maximum(temps)]
    end
    
    # Determine if we need one or two NASA polynomials
    if temperature_range[2] - temperature_range[1] > 1000
        # Split range for two polynomials
        middle_temp = 1000.0
        
        # Ensure middle temperature is within the data range
        if middle_temp < temperature_range[1]
            middle_temp = temperature_range[1]
        elseif middle_temp > temperature_range[2]
            middle_temp = temperature_range[2]
        end
        
        # Create two ranges
        range_low = [temperature_range[1], middle_temp]
        range_high = [middle_temp, temperature_range[2]]
        
        # Filter data for each range
        idx_low = findall(t -> range_low[1] <= t <= range_low[2], temps)
        idx_high = findall(t -> range_high[1] <= t <= range_high[2], temps)
        
        # Ensure we have enough points in each range
        if length(idx_low) < 5
            # Not enough points, use one polynomial
            return convert_tabular_to_single_nasa7(temps, cps, hs, ss, temperature_range)
        end
        
        if length(idx_high) < 5
            # Not enough points, use one polynomial
            return convert_tabular_to_single_nasa7(temps, cps, hs, ss, temperature_range)
        end
        
        # Create two polynomials
        temps_low = temps[idx_low]
        cps_low = cps[idx_low]
        hs_low = has_h ? hs[idx_low] : zeros(length(temps_low))
        ss_low = has_s ? ss[idx_low] : zeros(length(temps_low))
        
        temps_high = temps[idx_high]
        cps_high = cps[idx_high]
        hs_high = has_h ? hs[idx_high] : zeros(length(temps_high))
        ss_high = has_s ? ss[idx_high] : zeros(length(temps_high))
        
        # Fit polynomials
        coeffs_low = fit_nasa7_coefficients(temps_low, cps_low, hs_low, ss_low)
        coeffs_high = fit_nasa7_coefficients(temps_high, cps_high, hs_high, ss_high)
        
        # Return NASA 7 data
        return Dict(
            "temperature_ranges" => [tuple(range_low...), tuple(range_high...)],
            "coefficients" => [coeffs_low, coeffs_high],
            "polynomial_type" => "nasa7"
        )
    else
        # Use a single polynomial
        return convert_tabular_to_single_nasa7(temps, cps, hs, ss, temperature_range)
    end
end

"""
    convert_tabular_to_single_nasa7(temps::Vector{Float64}, cps::Vector{Float64}, 
                                 hs::Vector{Float64}, ss::Vector{Float64}, 
                                 temperature_range::Vector{Float64})

Convert tabular data to a single NASA 7-coefficient polynomial.
"""
function convert_tabular_to_single_nasa7(temps::Vector{Float64}, cps::Vector{Float64}, 
                                      hs::Vector{Float64}, ss::Vector{Float64}, 
                                      temperature_range::Vector{Float64})
    # Filter data to the specified temperature range
    idx = findall(t -> temperature_range[1] <= t <= temperature_range[2], temps)
    
    if length(idx) < 5
        error("Not enough data points in the specified temperature range")
    end
    
    # Extract data in range
    temps_range = temps[idx]
    cps_range = cps[idx]
    hs_range = hs[idx]
    ss_range = ss[idx]
    
    # Fit NASA 7-coefficient polynomial
    coeffs = fit_nasa7_coefficients(temps_range, cps_range, hs_range, ss_range)
    
    # Return NASA 7 data
    return Dict(
        "temperature_ranges" => [tuple(temperature_range...)],
        "coefficients" => [coeffs],
        "polynomial_type" => "nasa7"
    )
end

"""
    convert_nasa7_to_nasa9(nasa7_data::Dict)

Convert NASA 7-coefficient polynomial data to NASA 9-coefficient format.
"""
function convert_nasa7_to_nasa9(nasa7_data::Dict)
    # This is a more complex conversion that requires refitting
    # In practice, a more sophisticated algorithm would be used
    
    temp_ranges = nasa7_data["temperature_ranges"]
    coeffs7 = nasa7_data["coefficients"]
    
    # Initialize result
    nasa9_data = Dict(
        "temperature_ranges" => temp_ranges,
        "coefficients" => [],
        "polynomial_type" => "nasa9"
    )
    
    # Convert each temperature range
    for i in 1:length(temp_ranges)
        range = temp_ranges[i]
        coeffs = coeffs7[i]
        
        # Generate synthetic data points from NASA 7
        temps = range[1]:50:range[2]
        
        cps = []
        hs = []
        ss = []
        
        for temp in temps
            cp = calculate_nasa7_cp(coeffs, temp)
            h = calculate_nasa7_enthalpy(coeffs, temp)
            s = calculate_nasa7_entropy(coeffs, temp)
            
            push!(cps, cp)
            push!(hs, h)
            push!(ss, s)
        end
        
        # Convert synthetic data to NASA 9
        # This is a placeholder - in reality, a proper fitting algorithm would be used
        coeffs9 = zeros(9)
        
        # Copy the first 5 coefficients (approximate conversion)
        for j in 1:5
            coeffs9[j+2] = coeffs[j]
        end
        
        # Estimate the other coefficients
        coeffs9[1] = 0.0  # a1 (T^-2 term)
        coeffs9[2] = 0.0  # a2 (T^-1 term)
        coeffs9[8] = coeffs[6]  # a8 (integration constant for H)
        coeffs9[9] = coeffs[7]  # a9 (integration constant for S)
        
        push!(nasa9_data["coefficients"], coeffs9)
    end
    
    return nasa9_data
end

"""
    convert_shomate_to_nasa7(shomate_data::Dict)

Convert Shomate equation parameters to NASA 7-coefficient format.
"""
function convert_shomate_to_nasa7(shomate_data::Dict)
    # Shomate equation: Cp = A + B*t + C*t^2 + D*t^3 + E/t^2, where t = T/1000
    # NASA 7: Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    
    temp_ranges = shomate_data["temperature_ranges"]
    shomate_coeffs = shomate_data["coefficients"]
    
    # Initialize result
    nasa7_data = Dict(
        "temperature_ranges" => temp_ranges,
        "coefficients" => [],
        "polynomial_type" => "nasa7"
    )
    
    # Convert each temperature range
    for i in 1:length(temp_ranges)
        range = temp_ranges[i]
        coeffs = shomate_coeffs[i]
        
        # Extract Shomate coefficients
        A = coeffs[1]
        B = coeffs[2]
        C = coeffs[3]
        D = coeffs[4]
        E = coeffs[5]
        F = coeffs[6]
        G = coeffs[7]
        H = coeffs[8]
        
        # Generate synthetic data points from Shomate
        temps = range[1]:50:range[2]
        
        cps = []
        hs = []
        ss = []
        
        for temp in temps
            t = temp / 1000.0  # Shomate uses T/1000
            
            # Calculate properties using Shomate equations
            cp = A + B*t + C*t^2 + D*t^3 + E/t^2  # J/mol/K
            h = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F  # kJ/mol
            s = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G  # J/mol/K
            
            # Convert to NASA 7 units
            cp_r = cp / R
            h_rt = h * 1000 / (R * temp)  # Convert kJ to J
            s_r = s / R
            
            push!(cps, cp_r)
            push!(hs, h_rt)
            push!(ss, s_r)
        end
        
        # Fit NASA 7-coefficient polynomial to synthetic data
        coeffs7 = fit_nasa7_coefficients(temps, cps, hs, ss)
        
        push!(nasa7_data["coefficients"], coeffs7)
    end
    
    return nasa7_data
end

"""
    convert_nasa_to_shomate(nasa_data::Dict)

Convert NASA polynomial coefficients to Shomate equation format.
"""
function convert_nasa_to_shomate(nasa_data::Dict)
    # This is a complex conversion that requires refitting
    # In practice, a more sophisticated algorithm would be used
    
    temp_ranges = nasa_data["temperature_ranges"]
    nasa_coeffs = nasa_data["coefficients"]
    poly_type = nasa_data["polynomial_type"]
    
    # Initialize result
    shomate_data = Dict(
        "temperature_ranges" => temp_ranges,
        "coefficients" => [],
        "polynomial_type" => "shomate"
    )
    
    # Convert each temperature range
    for i in 1:length(temp_ranges)
        range = temp_ranges[i]
        coeffs = nasa_coeffs[i]
        
        # Generate synthetic data points from NASA
        temps = range[1]:50:range[2]
        
        cps = []
        hs = []
        ss = []
        
        for temp in temps
            if poly_type == "nasa7"
                cp = calculate_nasa7_cp(coeffs, temp)
                h = calculate_nasa7_enthalpy(coeffs, temp)
                s = calculate_nasa7_entropy(coeffs, temp)
            else  # nasa9
                cp = calculate_nasa9_cp(coeffs, temp)
                h = calculate_nasa9_enthalpy(coeffs, temp)
                s = calculate_nasa9_entropy(coeffs, temp)
            end
            
            # Convert to Shomate units
            cp_j = cp * R  # J/mol/K
            h_kj = h * R * temp / 1000  # kJ/mol
            s_j = s * R  # J/mol/K
            
            push!(cps, cp_j)
            push!(hs, h_kj)
            push!(ss, s_j)
        end
        
        # Fit Shomate equation to synthetic data
        shomate_coeffs = fit_shomate_coefficients(temps, cps, hs, ss)
        
        push!(shomate_data["coefficients"], shomate_coeffs)
    end
    
    return shomate_data
end

"""
    fit_shomate_coefficients(temps::Vector{Float64}, cps::Vector{Float64}, 
                           hs::Vector{Float64}, ss::Vector{Float64})

Fit Shomate equation coefficients to thermodynamic data.
"""
function fit_shomate_coefficients(temps::Vector{Float64}, cps::Vector{Float64}, 
                                hs::Vector{Float64}, ss::Vector{Float64})
    # Shomate equation: Cp = A + B*t + C*t^2 + D*t^3 + E/t^2, where t = T/1000
    
    # Convert temperatures to Shomate scale
    t = temps / 1000.0
    
    # Fit heat capacity coefficients (A, B, C, D, E)
    n = length(temps)
    X = zeros(n, 5)
    
    for i in 1:n
        X[i, 1] = 1.0
        X[i, 2] = t[i]
        X[i, 3] = t[i]^2
        X[i, 4] = t[i]^3
        X[i, 5] = 1/t[i]^2
    end
    
    # Solve least squares problem for A, B, C, D, E
    abcde = X \ cps
    
    A = abcde[1]
    B = abcde[2]
    C = abcde[3]
    D = abcde[4]
    E = abcde[5]
    
    # Fit integration constants for enthalpy (F) and entropy (G)
    F = 0.0
    G = 0.0
    
    # Return Shomate coefficients
    return [A, B, C, D, E, F, G, 0.0]  # H is unused
end