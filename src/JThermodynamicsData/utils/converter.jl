"""
    JThermodynamicsData.utils.converter

This module handles conversion between different thermodynamic polynomial formats.
- Converts between NASA-7 and NASA-9 polynomial formats
- Generates polynomial fits for thermodynamic data
- Handles uncertainty propagation in polynomial conversions
"""

using LsqFit
using Statistics
using DataFrames
using Plots
using Dates
using JSON

"""
    generate_nasa7_polynomial(temps::Vector{Float64}, cp_values::Vector{Float64}, 
                            h_values::Vector{Float64}, s_values::Vector{Float64},
                            temp_ranges::Vector{Vector{Float64}}=[[200.0, 1000.0], [1000.0, 6000.0]])

Generate NASA-7 polynomial coefficients from thermodynamic data.
"""
function generate_nasa7_polynomial(temps::Vector{Float64}, cp_values::Vector{Float64}, 
                                 h_values::Vector{Float64}, s_values::Vector{Float64},
                                 temp_ranges::Vector{Vector{Float64}}=[[200.0, 1000.0], [1000.0, 6000.0]])
    # Convert all values to dimensionless units (divide by R)
    cp_dimensionless = cp_values ./ 8.31446
    h_dimensionless = h_values ./ (8.31446 .* temps)  # H/RT
    s_dimensionless = s_values ./ 8.31446  # S/R
    
    coeffs = []
    
    for range_idx in 1:length(temp_ranges)
        range_min = temp_ranges[range_idx][1]
        range_max = temp_ranges[range_idx][2]
        
        # Filter data points in this range
        range_indices = findall(t -> range_min <= t <= range_max, temps)
        
        if isempty(range_indices)
            # No data in this range, use extrapolation or defaults
            push!(coeffs, [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            continue
        end
        
        range_temps = temps[range_indices]
        range_cp = cp_dimensionless[range_indices]
        range_h = h_dimensionless[range_indices]
        range_s = s_dimensionless[range_indices]
        
        # Initial guesses for coefficients
        initial_guess = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        # Fit Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        function cp_model(t, p)
            return p[1] .+ p[2] .* t .+ p[3] .* t.^2 .+ p[4] .* t.^3 .+ p[5] .* t.^4
        end
        
        # Using normalized temperatures for numerical stability
        norm_temps = range_temps ./ 1000.0
        
        # Fit Cp model
        try
            cp_fit = curve_fit(cp_model, norm_temps, range_cp, initial_guess[1:5])
            a1, a2, a3, a4, a5 = cp_fit.param
            
            # Use Cp coefficients to compute H and S coefficients
            # H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
            # S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
            
            # Compute a6 from enthalpy
            h_model(t, p) = cp_model(t, p[1:5]) .+ p[6] ./ t
            h_fit = curve_fit(h_model, norm_temps, range_h, [a1, a2, a3, a4, a5, 0.0])
            a6 = h_fit.param[6]
            
            # Compute a7 from entropy
            function s_model(t, p)
                return p[1] .* log.(t) .+ p[2] .* t .+ p[3] .* t.^2 ./ 2 .+ 
                       p[4] .* t.^3 ./ 3 .+ p[5] .* t.^4 ./ 4 .+ p[7]
            end
            s_fit = curve_fit(s_model, norm_temps, range_s, [a1, a2, a3, a4, a5, a6, 0.0])
            a7 = s_fit.param[7]
            
            # Rescale coefficients for unnormalized temperatures
            a2 = a2 / 1000.0
            a3 = a3 / 1000.0^2
            a4 = a4 / 1000.0^3
            a5 = a5 / 1000.0^4
            a6 = a6 * 1000.0
            
            push!(coeffs, [a1, a2, a3, a4, a5, a6, a7])
        catch e
            @warn "Error fitting NASA-7 polynomial: $e"
            push!(coeffs, [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        end
    end
    
    return coeffs
end

"""
    generate_nasa9_polynomial(temps::Vector{Float64}, cp_values::Vector{Float64}, 
                            h_values::Vector{Float64}, s_values::Vector{Float64},
                            temp_ranges::Vector{Vector{Float64}}=[[200.0, 1000.0], [1000.0, 6000.0], [6000.0, 20000.0]])

Generate NASA-9 polynomial coefficients from thermodynamic data.
"""
function generate_nasa9_polynomial(temps::Vector{Float64}, cp_values::Vector{Float64}, 
                                 h_values::Vector{Float64}, s_values::Vector{Float64},
                                 temp_ranges::Vector{Vector{Float64}}=[[200.0, 1000.0], [1000.0, 6000.0], [6000.0, 20000.0]])
    # Convert all values to dimensionless units (divide by R)
    cp_dimensionless = cp_values ./ 8.31446
    h_dimensionless = h_values ./ (8.31446 .* temps)  # H/RT
    s_dimensionless = s_values ./ 8.31446  # S/R
    
    coeffs = []
    
    for range_idx in 1:length(temp_ranges)
        range_min = temp_ranges[range_idx][1]
        range_max = temp_ranges[range_idx][2]
        
        # Filter data points in this range
        range_indices = findall(t -> range_min <= t <= range_max, temps)
        
        if isempty(range_indices)
            # No data in this range, use extrapolation or defaults
            push!(coeffs, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            continue
        end
        
        range_temps = temps[range_indices]
        range_cp = cp_dimensionless[range_indices]
        range_h = h_dimensionless[range_indices]
        range_s = s_dimensionless[range_indices]
        
        # Initial guesses for coefficients
        initial_guess = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        # Fit Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4 + a6*T^5 + a7*T^6 + a8*T^7
        function cp_model(t, p)
            return p[1] .+ p[2] .* t .+ p[3] .* t.^2 .+ p[4] .* t.^3 .+ 
                   p[5] .* t.^4 .+ p[6] .* t.^5 .+ p[7] .* t.^6 .+ p[8] .* t.^7
        end
        
        # Using normalized temperatures for numerical stability
        norm_temps = range_temps ./ 1000.0
        
        # Fit Cp model
        try
            cp_fit = curve_fit(cp_model, norm_temps, range_cp, initial_guess[1:8])
            a1, a2, a3, a4, a5, a6, a7, a8 = cp_fit.param
            
            # Use Cp coefficients to compute H and S coefficients
            # a9 (integration constant) for entropy
            
            # Compute a9 from entropy
            function s_model(t, p)
                return p[1] .* log.(t) .+ p[2] .* t .+ p[3] .* t.^2 ./ 2 .+ 
                       p[4] .* t.^3 ./ 3 .+ p[5] .* t.^4 ./ 4 .+ 
                       p[6] .* t.^5 ./ 5 .+ p[7] .* t.^6 ./ 6 .+ 
                       p[8] .* t.^7 ./ 7 .+ p[9]
            end
            s_fit = curve_fit(s_model, norm_temps, range_s, [a1, a2, a3, a4, a5, a6, a7, a8, 0.0])
            a9 = s_fit.param[9]
            
            # Rescale coefficients for unnormalized temperatures
            a2 = a2 / 1000.0
            a3 = a3 / 1000.0^2
            a4 = a4 / 1000.0^3
            a5 = a5 / 1000.0^4
            a6 = a6 / 1000.0^5
            a7 = a7 / 1000.0^6
            a8 = a8 / 1000.0^7
            
            push!(coeffs, [a1, a2, a3, a4, a5, a6, a7, a8, a9])
        catch e
            @warn "Error fitting NASA-9 polynomial: $e"
            push!(coeffs, [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        end
    end
    
    return coeffs
end

"""
    nasa7_to_nasa9(nasa7_coeffs::Vector{Vector{Float64}})

Convert NASA-7 polynomial coefficients to NASA-9 format.
"""
function nasa7_to_nasa9(nasa7_coeffs::Vector{Vector{Float64}})
    nasa9_coeffs = []
    
    for coeffs in nasa7_coeffs
        a1, a2, a3, a4, a5, a6, a7 = coeffs
        
        # NASA-9 uses 9 coefficients, with a6 and a7 for high-temperature terms
        # We'll set these to 0 for simple conversion
        push!(nasa9_coeffs, [a1, a2, a3, a4, a5, 0.0, 0.0, 0.0, a7])
    end
    
    return nasa9_coeffs
end

"""
    nasa9_to_nasa7(nasa9_coeffs::Vector{Vector{Float64}})

Convert NASA-9 polynomial coefficients to NASA-7 format.
This is lossy - the high order terms are dropped.
"""
function nasa9_to_nasa7(nasa9_coeffs::Vector{Vector{Float64}})
    nasa7_coeffs = []
    
    for coeffs in nasa9_coeffs
        a1, a2, a3, a4, a5 = coeffs[1:5]
        a9 = coeffs[9]  # Integration constant for entropy
        
        # For h, we need to introduce a6 term (b2 in NASA-7)
        a6 = 0.0  # Estimate this from enthalpy data if available
        
        push!(nasa7_coeffs, [a1, a2, a3, a4, a5, a6, a9])
    end
    
    return nasa7_coeffs
end

"""
    generate_polynomial_with_uncertainty(temps::Vector{Float64}, property_values::Dict, 
                                      temp_ranges::Vector{Vector{Float64}})

Generate both NASA-7 and NASA-9 polynomials with uncertainty estimates.
"""
function generate_polynomial_with_uncertainty(temps::Vector{Float64}, property_values::Dict, 
                                           temp_ranges::Vector{Vector{Float64}})
    # Extract data
    cp_values = property_values["Cp"]["values"]
    h_values = property_values["H"]["values"]
    s_values = property_values["S"]["values"]
    
    cp_uncertainties = property_values["Cp"]["uncertainties"]
    h_uncertainties = property_values["H"]["uncertainties"]
    s_uncertainties = property_values["S"]["uncertainties"]
    
    # Generate polynomials for mean values
    nasa7_mean = generate_nasa7_polynomial(temps, cp_values, h_values, s_values, temp_ranges)
    nasa9_mean = generate_nasa9_polynomial(temps, cp_values, h_values, s_values, temp_ranges)
    
    # Generate polynomials for upper bound (mean + uncertainty)
    nasa7_upper = generate_nasa7_polynomial(
        temps, 
        cp_values .+ cp_uncertainties, 
        h_values .+ h_uncertainties, 
        s_values .+ s_uncertainties, 
        temp_ranges
    )
    
    nasa9_upper = generate_nasa9_polynomial(
        temps, 
        cp_values .+ cp_uncertainties, 
        h_values .+ h_uncertainties, 
        s_values .+ s_uncertainties, 
        temp_ranges
    )
    
    # Generate polynomials for lower bound (mean - uncertainty)
    nasa7_lower = generate_nasa7_polynomial(
        temps, 
        cp_values .- cp_uncertainties, 
        h_values .- h_uncertainties, 
        s_values .- s_uncertainties, 
        temp_ranges
    )
    
    nasa9_lower = generate_nasa9_polynomial(
        temps, 
        cp_values .- cp_uncertainties, 
        h_values .- h_uncertainties, 
        s_values .- s_uncertainties, 
        temp_ranges
    )
    
    # Calculate uncertainty polynomials (half the difference between upper and lower)
    nasa7_uncertainty = []
    nasa9_uncertainty = []
    
    for i in 1:length(nasa7_mean)
        # For each coefficient, uncertainty is half the difference between upper and lower
        nasa7_unc = (nasa7_upper[i] .- nasa7_lower[i]) ./ 2
        nasa9_unc = (nasa9_upper[i] .- nasa9_lower[i]) ./ 2
        
        push!(nasa7_uncertainty, nasa7_unc)
        push!(nasa9_uncertainty, nasa9_unc)
    end
    
    return Dict(
        "nasa7" => Dict(
            "mean" => nasa7_mean,
            "uncertainty" => nasa7_uncertainty
        ),
        "nasa9" => Dict(
            "mean" => nasa9_mean,
            "uncertainty" => nasa9_uncertainty
        )
    )
end

"""
    export_polynomial_data(species_name::String, formula::String, cas::String, polynomials::Dict, 
                        temp_ranges::Vector{Vector{Float64}}, output_dir::String)

Export polynomial data in various formats (JSON, CSV, formatted text).
"""
function export_polynomial_data(species_name::String, formula::String, cas::String, polynomials::Dict, 
                             temp_ranges::Vector{Vector{Float64}}, output_dir::String)
    # Create output directory if needed
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Create subdirectories
    json_dir = joinpath(output_dir, "json")
    csv_dir = joinpath(output_dir, "csv")
    text_dir = joinpath(output_dir, "text")
    
    for dir in [json_dir, csv_dir, text_dir]
        if !isdir(dir)
            mkpath(dir)
        end
    end
    
    # Export as JSON
    json_data = Dict(
        "species_name" => species_name,
        "formula" => formula,
        "cas_registry" => cas,
        "temperature_ranges" => temp_ranges,
        "nasa7" => polynomials["nasa7"],
        "nasa9" => polynomials["nasa9"],
        "metadata" => Dict(
            "generated" => string(now()),
            "version" => "1.0.0"
        )
    )
    
    open(joinpath(json_dir, "$(species_name)_polynomials.json"), "w") do io
        JSON.print(io, json_data, 4)  # Pretty print with 4-space indent
    end
    
    # Export as CSV - Create table with temperature, property values, and uncertainties
    temps = collect(temp_ranges[1][1]:10:temp_ranges[end][2])
    
    # Calculate properties from polynomials
    cp_values = []
    h_values = []
    s_values = []
    g_values = []
    
    cp_uncertainties = []
    h_uncertainties = []
    s_uncertainties = []
    g_uncertainties = []
    
    for t in temps
        # Find applicable range
        range_idx = 1
        for (i, range) in enumerate(temp_ranges)
            if range[1] <= t <= range[2]
                range_idx = i
                break
            end
        end
        
        # Get coefficients
        nasa7_coeffs = polynomials["nasa7"]["mean"][range_idx]
        nasa7_unc = polynomials["nasa7"]["uncertainty"][range_idx]
        
        # Calculate properties using NASA-7 formula
        cp = calculate_nasa7_cp(nasa7_coeffs, t)
        h = calculate_nasa7_enthalpy(nasa7_coeffs, t)
        s = calculate_nasa7_entropy(nasa7_coeffs, t)
        g = h - t * s / 1000  # kJ/mol
        
        # Calculate uncertainties
        cp_unc = calculate_nasa7_cp_uncertainty(nasa7_coeffs, nasa7_unc, t)
        h_unc = calculate_nasa7_enthalpy_uncertainty(nasa7_coeffs, nasa7_unc, t)
        s_unc = calculate_nasa7_entropy_uncertainty(nasa7_coeffs, nasa7_unc, t)
        g_unc = sqrt(h_unc^2 + (t * s_unc / 1000)^2)  # Propagate uncertainty
        
        push!(cp_values, cp)
        push!(h_values, h)
        push!(s_values, s)
        push!(g_values, g)
        
        push!(cp_uncertainties, cp_unc)
        push!(h_uncertainties, h_unc)
        push!(s_uncertainties, s_unc)
        push!(g_uncertainties, g_unc)
    end
    
    # Create DataFrame
    df = DataFrame(
        Temperature_K = temps,
        Cp_JmolK = cp_values,
        Cp_Uncertainty_JmolK = cp_uncertainties,
        H_kJmol = h_values,
        H_Uncertainty_kJmol = h_uncertainties,
        S_JmolK = s_values,
        S_Uncertainty_JmolK = s_uncertainties,
        G_kJmol = g_values,
        G_Uncertainty_kJmol = g_uncertainties
    )
    
    # Save CSV
    CSV.write(joinpath(csv_dir, "$(species_name)_thermodynamic_data.csv"), df)
    
    # Export as formatted text
    open(joinpath(text_dir, "$(species_name)_nasa7.txt"), "w") do io
        println(io, "")
        println(io, "NASA-7 Polynomial Coefficients for $(species_name) ($(formula)):")
        
        # Format temperature ranges
        ranges_str = join(["$(range[1])-$(range[2]) K" for range in temp_ranges], " and ")
        println(io, "Temperature ranges: $(ranges_str)")
        println(io, "")
        println(io, "! NASA-7 polynomial in standard format")
        
        # Format 7-term NASA polynomials
        nasa7_coeffs = polynomials["nasa7"]["mean"]
        
        # First line with species info
        println(io, "$(rpad(species_name, 16))$(rpad(formula, 16))$(lpad(string(Int(temp_ranges[1][1])), 8)).00$(lpad(string(Int(temp_ranges[end][2])), 9)).00$(lpad(string(Int(temp_ranges[1][2])), 9)).00")
        
        # High-temperature coefficients
        high_temp_coeffs = nasa7_coeffs[end]
        println(io, "$(lpad(@sprintf("%.8e", high_temp_coeffs[1]), 16))$(lpad(@sprintf("%.8e", high_temp_coeffs[2]), 16))$(lpad(@sprintf("%.8e", high_temp_coeffs[3]), 16))$(lpad(@sprintf("%.8e", high_temp_coeffs[4]), 16))$(lpad(@sprintf("%.8e", high_temp_coeffs[5]), 16))    2")
        println(io, "$(lpad(@sprintf("%.8e", high_temp_coeffs[6]), 16))$(lpad(@sprintf("%.8e", high_temp_coeffs[7]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][1]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][2]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][3]), 16))    3")
        println(io, "$(lpad(@sprintf("%.8e", nasa7_coeffs[1][4]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][5]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][6]), 16))$(lpad(@sprintf("%.8e", nasa7_coeffs[1][7]), 16))                   4")
    end
    
    open(joinpath(text_dir, "$(species_name)_nasa9.txt"), "w") do io
        println(io, "")
        println(io, "NASA-9 Polynomial Coefficients for $(species_name) ($(formula)):")
        
        # Format temperature ranges
        ranges_str = join(["$(range[1])-$(range[2]) K" for range in temp_ranges], ", ")
        println(io, "Temperature ranges: $(ranges_str)")
        println(io, "")
        println(io, "! NASA-9 polynomial in standard format")
        
        # Format 9-term NASA polynomials
        nasa9_coeffs = polynomials["nasa9"]["mean"]
        
        # Format NASA-9 polynomial
        # Format varies by implementation, this is a basic version
        println(io, "$(species_name) $(formula)")
        println(io, "! CAS: $(cas)")
        println(io, "! Temperature ranges: $(ranges_str)")
        
        for (i, range) in enumerate(temp_ranges)
            coeffs = nasa9_coeffs[i]
            println(io, "! Range $(i): $(range[1])-$(range[2]) K")
            println(io, "$(lpad(@sprintf("%.8e", coeffs[1]), 16))$(lpad(@sprintf("%.8e", coeffs[2]), 16))$(lpad(@sprintf("%.8e", coeffs[3]), 16))$(lpad(@sprintf("%.8e", coeffs[4]), 16))")
            println(io, "$(lpad(@sprintf("%.8e", coeffs[5]), 16))$(lpad(@sprintf("%.8e", coeffs[6]), 16))$(lpad(@sprintf("%.8e", coeffs[7]), 16))$(lpad(@sprintf("%.8e", coeffs[8]), 16))")
            println(io, "$(lpad(@sprintf("%.8e", coeffs[9]), 16))")
        end
    end
    
    return Dict(
        "json_path" => joinpath(json_dir, "$(species_name)_polynomials.json"),
        "csv_path" => joinpath(csv_dir, "$(species_name)_thermodynamic_data.csv"),
        "nasa7_path" => joinpath(text_dir, "$(species_name)_nasa7.txt"),
        "nasa9_path" => joinpath(text_dir, "$(species_name)_nasa9.txt")
    )
end

# Helper functions for calculating thermodynamic properties from NASA-7 polynomials

"""
    calculate_nasa7_cp(coeffs::Vector{Float64}, temperature::Float64)

Calculate heat capacity (Cp) from NASA-7 polynomial coefficients.
"""
function calculate_nasa7_cp(coeffs::Vector{Float64}, temperature::Float64)
    a1, a2, a3, a4, a5 = coeffs[1:5]
    
    cp_r = a1 + a2*temperature + a3*temperature^2 + a4*temperature^3 + a5*temperature^4
    
    # Convert to J/mol/K (R = 8.31446 J/mol/K)
    return cp_r * 8.31446
end

"""
    calculate_nasa7_enthalpy(coeffs::Vector{Float64}, temperature::Float64)

Calculate enthalpy (H) from NASA-7 polynomial coefficients.
"""
function calculate_nasa7_enthalpy(coeffs::Vector{Float64}, temperature::Float64)
    a1, a2, a3, a4, a5, a6 = coeffs[1:6]
    
    h_rt = a1 + a2*temperature/2 + a3*temperature^2/3 + a4*temperature^3/4 + a5*temperature^4/5 + a6/temperature
    
    # Convert to kJ/mol (R*T / 1000)
    return h_rt * 8.31446 * temperature / 1000
end

"""
    calculate_nasa7_entropy(coeffs::Vector{Float64}, temperature::Float64)

Calculate entropy (S) from NASA-7 polynomial coefficients.
"""
function calculate_nasa7_entropy(coeffs::Vector{Float64}, temperature::Float64)
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    
    s_r = a1*log(temperature) + a2*temperature + a3*temperature^2/2 + a4*temperature^3/3 + a5*temperature^4/4 + a7
    
    # Convert to J/mol/K
    return s_r * 8.31446
end

"""
    calculate_nasa7_cp_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)

Calculate uncertainty in heat capacity from NASA-7 polynomial coefficient uncertainties.
"""
function calculate_nasa7_cp_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)
    # Use error propagation formula
    variance = unc[1]^2 +
               (unc[2] * temperature)^2 +
               (unc[3] * temperature^2)^2 +
               (unc[4] * temperature^3)^2 +
               (unc[5] * temperature^4)^2
    
    # Convert to J/mol/K
    return sqrt(variance) * 8.31446
end

"""
    calculate_nasa7_enthalpy_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)

Calculate uncertainty in enthalpy from NASA-7 polynomial coefficient uncertainties.
"""
function calculate_nasa7_enthalpy_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)
    # Use error propagation formula
    variance = unc[1]^2 +
               (unc[2] * temperature/2)^2 +
               (unc[3] * temperature^2/3)^2 +
               (unc[4] * temperature^3/4)^2 +
               (unc[5] * temperature^4/5)^2 +
               (unc[6] / temperature)^2
    
    # Convert to kJ/mol
    return sqrt(variance) * 8.31446 * temperature / 1000
end

"""
    calculate_nasa7_entropy_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)

Calculate uncertainty in entropy from NASA-7 polynomial coefficient uncertainties.
"""
function calculate_nasa7_entropy_uncertainty(coeffs::Vector{Float64}, unc::Vector{Float64}, temperature::Float64)
    # Use error propagation formula
    variance = (unc[1] * log(temperature))^2 +
               (unc[2] * temperature)^2 +
               (unc[3] * temperature^2/2)^2 +
               (unc[4] * temperature^3/3)^2 +
               (unc[5] * temperature^4/4)^2 +
               unc[7]^2
    
    # Convert to J/mol/K
    return sqrt(variance) * 8.31446
end