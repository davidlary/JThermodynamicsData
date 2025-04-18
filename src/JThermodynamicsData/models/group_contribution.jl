"""
Group contribution methods for estimating thermodynamic properties.
"""

"""
    estimate_properties_group_contribution(formula::String, temperature::Float64)

Estimate thermodynamic properties using group contribution methods.
This function selects the most appropriate group contribution method based on the molecule.
"""
function estimate_properties_group_contribution(formula::String, temperature::Float64)
    # Try Benson group additivity first (more accurate)
    try
        return benson_group_additivity(formula, temperature)
    catch e
        # Fall back to simpler group contribution method
        @warn "Benson group additivity failed for $formula: $e. Using simpler method."
        return simple_group_contribution(formula, temperature)
    end
end

"""
    simple_group_contribution(formula::String, temperature::Float64)

A simplified group contribution method based on element counts.
"""
function simple_group_contribution(formula::String, temperature::Float64)
    # Parse formula to identify elements
    elements = extract_formula_elements(formula)
    
    if length(elements) == 0
        error("Could not parse formula: $formula")
    end
    
    # Simple group contribution method based on element counts
    
    # Heat capacity contribution (J/mol/K)
    cp_contrib = Dict(
        "C" => 8.0,   # Carbon
        "H" => 2.0,   # Hydrogen
        "O" => 6.0,   # Oxygen
        "N" => 6.0,   # Nitrogen
        "S" => 6.0,   # Sulfur
        "Cl" => 8.0,  # Chlorine
        "F" => 5.0,   # Fluorine
        "Br" => 8.5,  # Bromine
        "I" => 9.0    # Iodine
    )
    
    # Entropy contribution (J/mol/K)
    s_contrib = Dict(
        "C" => 40.0,  # Carbon
        "H" => 7.0,   # Hydrogen
        "O" => 25.0,  # Oxygen
        "N" => 20.0,  # Nitrogen
        "S" => 30.0,  # Sulfur
        "Cl" => 35.0, # Chlorine
        "F" => 30.0,  # Fluorine
        "Br" => 40.0, # Bromine
        "I" => 42.0   # Iodine
    )
    
    # Enthalpy of formation contribution (kJ/mol)
    h_contrib = Dict(
        "C" => -10.0,  # Carbon
        "H" => -4.0,   # Hydrogen
        "O" => -20.0,  # Oxygen
        "N" => 10.0,   # Nitrogen
        "S" => 15.0,   # Sulfur
        "Cl" => -5.0,  # Chlorine
        "F" => -25.0,  # Fluorine
        "Br" => 0.0,   # Bromine
        "I" => 10.0    # Iodine
    )
    
    # Calculate base properties at 298.15 K
    cp_298 = 0.0
    s_298 = 0.0
    h_298 = 0.0
    
    for (element, count) in elements
        if haskey(cp_contrib, element)
            cp_298 += cp_contrib[element] * count
        end
        
        if haskey(s_contrib, element)
            s_298 += s_contrib[element] * count
        end
        
        if haskey(h_contrib, element)
            h_298 += h_contrib[element] * count
        end
    end
    
    # Adjust for temperature dependence (simplified model)
    # Cp = a + b*T + c*T^2
    # Where a, b, c are estimated from element counts
    a = cp_298
    b = sum([0.1 * count for (_, count) in elements])
    c = -sum([0.0001 * count for (_, count) in elements])
    
    # Calculate Cp at the requested temperature
    cp = a + b * (temperature - 298.15) + c * (temperature - 298.15)^2
    
    # Calculate enthalpy change from 298.15 K to T
    h_change = a * (temperature - 298.15) + b * (temperature - 298.15)^2 / 2 + c * (temperature - 298.15)^3 / 3
    h = h_298 + h_change
    
    # Calculate entropy change from 298.15 K to T
    s_change = a * log(temperature / 298.15) + b * (temperature - 298.15) + c * (temperature - 298.15)^2 / 2
    s = s_298 + s_change
    
    # Calculate Gibbs free energy
    g = h - temperature * s / 1000  # Convert S from J/mol/K to kJ/mol/K
    
    # Estimate uncertainties (rough estimates)
    cp_uncertainty = max(cp * 0.2, 5.0)  # 20% or at least 5 J/mol/K
    h_uncertainty = max(abs(h) * 0.3, 20.0)  # 30% or at least 20 kJ/mol
    s_uncertainty = max(s * 0.25, 10.0)  # 25% or at least 10 J/mol/K
    g_uncertainty = max(abs(g) * 0.3, 20.0)  # 30% or at least 20 kJ/mol
    
    # Create result with uncertainties
    return Dict(
        "temperature" => temperature,
        "Cp" => cp,
        "Cp_uncertainty" => cp_uncertainty,
        "H" => h,
        "H_uncertainty" => h_uncertainty,
        "S" => s,
        "S_uncertainty" => s_uncertainty,
        "G" => g,
        "G_uncertainty" => g_uncertainty,
        "method" => "simple_group_contribution"
    )
end

"""
    benson_group_additivity(formula::String, temperature::Float64)

Estimate thermodynamic properties using Benson's group additivity method.
This is one of the most accurate group contribution methods for gas-phase thermochemistry.
"""
function benson_group_additivity(formula::String, temperature::Float64)
    # Parse formula to identify elements
    elements = extract_formula_elements(formula)
    
    if length(elements) == 0
        error("Could not parse formula: $formula")
    end
    
    # Calculate molecular weight
    mw = 0.0
    for (element, count) in elements
        if haskey(ATOMIC_MASSES, element)
            mw += ATOMIC_MASSES[element] * count
        else
            @warn "Unknown element in formula: $element"
        end
    end
    
    # Determine molecule type based on elements
    is_hydrocarbon = haskey(elements, "C") && haskey(elements, "H") && length(elements) <= 2
    has_oxygen = haskey(elements, "O")
    has_nitrogen = haskey(elements, "N")
    has_halogen = any(haskey(elements, x) for x in ["F", "Cl", "Br", "I"])
    
    # Benson group contributions for common groups
    # Values for [Cp(298K), H(298K), S(298K)]
    # Units: Cp in J/mol/K, H in kJ/mol, S in J/mol/K
    group_contrib = Dict(
        # Alkane groups
        "C-(C)(H)3" => [35.1, -42.9, 127.2],  # Primary carbon (methyl)
        "C-(C)2(H)2" => [26.2, -20.9, 79.5],  # Secondary carbon
        "C-(C)3(H)" => [22.5, -8.8, 40.3],    # Tertiary carbon
        "C-(C)4" => [15.2, 4.5, -33.6],       # Quaternary carbon
        
        # Alkene groups
        "C=(C)(H)2" => [21.8, 26.9, 115.9],   # =CH2
        "C=(C)2(H)" => [19.6, 35.0, 40.0],    # =CH-
        "C=(C)3" => [18.0, 50.0, -32.0],      # =C<
        
        # Alkyne groups
        "C≡(C)(H)" => [17.4, 115.0, 25.0],    # ≡CH
        "C≡(C)2" => [23.5, 113.0, 8.0],       # ≡C-
        
        # Oxygen-containing groups
        "O-(C)(H)" => [18.8, -186.0, 112.0],  # Alcohol OH
        "O-(C)2" => [13.9, -120.0, 22.0],     # Ether O
        "C=O-(C)(H)" => [34.8, -146.0, 167.0], # Aldehyde
        "C=O-(C)2" => [31.0, -133.0, 150.0],  # Ketone
        "C=O-(O)(C)" => [31.0, -326.0, 169.0], # Ester
        "C=O-(O)(H)" => [35.1, -387.0, 290.0], # Carboxylic acid
        
        # Nitrogen-containing groups
        "N-(C)(H)2" => [26.5, 32.8, 196.0],   # Primary amine
        "N-(C)2(H)" => [23.0, 55.0, 114.0],   # Secondary amine
        "N-(C)3" => [19.0, 88.0, 39.0],       # Tertiary amine
        "N=(C)(C)" => [22.5, 68.0, 112.0],    # Imine
        "C≡N" => [23.0, 164.0, 94.0],         # Nitrile
        
        # Halogen-containing groups
        "C-(F)(C)(H)2" => [25.5, -220.0, 140.0], # CH2F
        "C-(Cl)(C)(H)2" => [31.0, -61.0, 150.0], # CH2Cl
        "C-(Br)(C)(H)2" => [31.0, -32.0, 155.0], # CH2Br
        "C-(I)(C)(H)2" => [31.5, 8.0, 160.0],    # CH2I
        
        # Ring corrections (additive)
        "ring-C3" => [27.0, 115.0, -155.0],   # Cyclopropane
        "ring-C4" => [27.0, 110.0, -120.0],   # Cyclobutane
        "ring-C5" => [7.0, 25.0, -108.0],     # Cyclopentane
        "ring-C6" => [-6.0, 0.0, -121.0],     # Cyclohexane
        
        # Aromatic groups
        "C-(H)(CB)2" => [15.8, -17.3, 42.4],  # Aromatic CH
        "C-(C)(CB)2" => [16.0, -9.0, -33.0],  # Substituted aromatic C
        
        # Default element contributions (for unrecognized groups)
        "C" => [16.0, 0.0, 40.0],             # Carbon
        "H" => [7.0, 0.0, 27.0],              # Hydrogen
        "O" => [7.5, -100.0, 45.0],           # Oxygen
        "N" => [8.0, 20.0, 50.0],             # Nitrogen
        "F" => [5.0, -150.0, 35.0],           # Fluorine
        "Cl" => [8.0, -35.0, 40.0],           # Chlorine
        "Br" => [8.5, -15.0, 48.0],           # Bromine
        "I" => [9.0, 10.0, 55.0]              # Iodine
    )
    
    # Temperature-dependent coefficients for Cp
    # Format: [a, b, c, d] for Cp = a + b*T + c*T^2 + d*T^3
    cp_coeffs = Dict(
        "C-(C)(H)3" => [3.5, 0.02, -1.0e-5, 1.5e-9],
        "C-(C)2(H)2" => [2.5, 0.015, -8.0e-6, 1.2e-9],
        "C-(C)3(H)" => [2.2, 0.01, -5.0e-6, 8.0e-10],
        "C-(C)4" => [1.5, 0.008, -4.0e-6, 6.0e-10],
        "C=(C)(H)2" => [2.0, 0.015, -7.0e-6, 1.0e-9],
        "C=(C)2(H)" => [1.8, 0.012, -6.0e-6, 9.0e-10],
        "C=(C)3" => [1.7, 0.01, -5.0e-6, 7.0e-10],
        "O-(C)(H)" => [1.9, 0.01, -5.0e-6, 7.0e-10],
        "O-(C)2" => [1.4, 0.008, -4.0e-6, 6.0e-10],
        "C=O-(C)(H)" => [3.4, 0.015, -7.0e-6, 1.0e-9],
        "C=O-(C)2" => [3.0, 0.012, -6.0e-6, 9.0e-10],
        "C=O-(O)(C)" => [3.0, 0.012, -6.0e-6, 9.0e-10],
        "C=O-(O)(H)" => [3.4, 0.015, -7.0e-6, 1.0e-9],
        "C" => [1.5, 0.01, -5.0e-6, 7.0e-10],
        "H" => [0.7, 0.001, -5.0e-7, 7.0e-11],
        "O" => [0.8, 0.003, -1.5e-6, 2.0e-10],
        "N" => [0.8, 0.005, -2.5e-6, 3.5e-10],
        "F" => [0.5, 0.001, -5.0e-7, 7.0e-11],
        "Cl" => [0.8, 0.001, -5.0e-7, 7.0e-11],
        "Br" => [0.85, 0.001, -5.0e-7, 7.0e-11],
        "I" => [0.9, 0.001, -5.0e-7, 7.0e-11]
    )
    
    # Use group contributions based on the approximate structure
    # This is a simplified approach - a real implementation would
    # require a proper structure or SMILES representation
    
    # Initialize properties
    cp_298 = 0.0
    h_298 = 0.0
    s_298 = 0.0
    
    # Initialize Cp coefficients
    a_sum = 0.0
    b_sum = 0.0
    c_sum = 0.0
    d_sum = 0.0
    
    # Process hydrocarbon (simplest case)
    if is_hydrocarbon
        c_count = elements["C"]
        h_count = elements["H"]
        
        # Identify the structure based on C:H ratio
        if c_count == 1 && h_count == 4
            # Methane (CH4)
            group = "C-(C)(H)3"  # Using methyl group for methane
            cp_298 += group_contrib[group][1]
            h_298 += group_contrib[group][2]
            s_298 += group_contrib[group][3]
            
            # Add Cp coefficients
            if haskey(cp_coeffs, group)
                a_sum += cp_coeffs[group][1]
                b_sum += cp_coeffs[group][2]
                c_sum += cp_coeffs[group][3]
                d_sum += cp_coeffs[group][4]
            end
        elseif c_count == 2 && h_count == 6
            # Ethane (C2H6)
            group = "C-(C)(H)3"  # Two methyl groups
            cp_298 += 2 * group_contrib[group][1]
            h_298 += 2 * group_contrib[group][2]
            s_298 += 2 * group_contrib[group][3]
            
            # Add Cp coefficients
            if haskey(cp_coeffs, group)
                a_sum += 2 * cp_coeffs[group][1]
                b_sum += 2 * cp_coeffs[group][2]
                c_sum += 2 * cp_coeffs[group][3]
                d_sum += 2 * cp_coeffs[group][4]
            end
        elseif c_count == 2 && h_count == 4
            # Ethylene (C2H4)
            group = "C=(C)(H)2"  # Two =CH2 groups
            cp_298 += 2 * group_contrib[group][1]
            h_298 += 2 * group_contrib[group][2]
            s_298 += 2 * group_contrib[group][3]
            
            # Add Cp coefficients
            if haskey(cp_coeffs, group)
                a_sum += 2 * cp_coeffs[group][1]
                b_sum += 2 * cp_coeffs[group][2]
                c_sum += 2 * cp_coeffs[group][3]
                d_sum += 2 * cp_coeffs[group][4]
            end
        elseif c_count == 2 && h_count == 2
            # Acetylene (C2H2)
            group = "C≡(C)(H)"  # Two ≡CH groups
            cp_298 += 2 * group_contrib[group][1]
            h_298 += 2 * group_contrib[group][2]
            s_298 += 2 * group_contrib[group][3]
            
            # Add Cp coefficients
            if haskey(cp_coeffs, group)
                a_sum += 2 * cp_coeffs[group][1]
                b_sum += 2 * cp_coeffs[group][2]
                c_sum += 2 * cp_coeffs[group][3]
                d_sum += 2 * cp_coeffs[group][4]
            end
        elseif c_count >= 3 && h_count >= 8  # Longer alkanes
            # Estimate as a mixture of methyl and methylene groups
            methyl_count = 2  # Two end groups
            methylene_count = c_count - 2  # Internal groups
            
            # Add methyl groups
            group = "C-(C)(H)3"
            cp_298 += methyl_count * group_contrib[group][1]
            h_298 += methyl_count * group_contrib[group][2]
            s_298 += methyl_count * group_contrib[group][3]
            
            # Add Cp coefficients for methyl
            if haskey(cp_coeffs, group)
                a_sum += methyl_count * cp_coeffs[group][1]
                b_sum += methyl_count * cp_coeffs[group][2]
                c_sum += methyl_count * cp_coeffs[group][3]
                d_sum += methyl_count * cp_coeffs[group][4]
            end
            
            # Add methylene groups
            group = "C-(C)2(H)2"
            cp_298 += methylene_count * group_contrib[group][1]
            h_298 += methylene_count * group_contrib[group][2]
            s_298 += methylene_count * group_contrib[group][3]
            
            # Add Cp coefficients for methylene
            if haskey(cp_coeffs, group)
                a_sum += methylene_count * cp_coeffs[group][1]
                b_sum += methylene_count * cp_coeffs[group][2]
                c_sum += methylene_count * cp_coeffs[group][3]
                d_sum += methylene_count * cp_coeffs[group][4]
            end
            
            # Check if it's likely cyclic
            if h_count <= 2 * c_count
                # Add ring correction based on size
                ring_type = "ring-C$c_count"
                if haskey(group_contrib, ring_type)
                    cp_298 += group_contrib[ring_type][1]
                    h_298 += group_contrib[ring_type][2]
                    s_298 += group_contrib[ring_type][3]
                else
                    # Use C6 ring as default for larger rings
                    cp_298 += group_contrib["ring-C6"][1]
                    h_298 += group_contrib["ring-C6"][2]
                    s_298 += group_contrib["ring-C6"][3]
                end
            end
        else
            # For more complex hydrocarbons, use element-based estimation
            cp_298 += elements["C"] * group_contrib["C"][1]
            h_298 += elements["C"] * group_contrib["C"][2]
            s_298 += elements["C"] * group_contrib["C"][3]
            
            cp_298 += elements["H"] * group_contrib["H"][1]
            h_298 += elements["H"] * group_contrib["H"][2]
            s_298 += elements["H"] * group_contrib["H"][3]
            
            # Add Cp coefficients
            a_sum += elements["C"] * cp_coeffs["C"][1]
            b_sum += elements["C"] * cp_coeffs["C"][2]
            c_sum += elements["C"] * cp_coeffs["C"][3]
            d_sum += elements["C"] * cp_coeffs["C"][4]
            
            a_sum += elements["H"] * cp_coeffs["H"][1]
            b_sum += elements["H"] * cp_coeffs["H"][2]
            c_sum += elements["H"] * cp_coeffs["H"][3]
            d_sum += elements["H"] * cp_coeffs["H"][4]
        end
    else
        # For non-hydrocarbons, use element-based estimation
        for (element, count) in elements
            if haskey(group_contrib, element)
                cp_298 += count * group_contrib[element][1]
                h_298 += count * group_contrib[element][2]
                s_298 += count * group_contrib[element][3]
                
                # Add Cp coefficients
                if haskey(cp_coeffs, element)
                    a_sum += count * cp_coeffs[element][1]
                    b_sum += count * cp_coeffs[element][2]
                    c_sum += count * cp_coeffs[element][3]
                    d_sum += count * cp_coeffs[element][4]
                end
            end
        end
        
        # Apply corrections for certain functional groups
        if has_oxygen
            # Check for common oxygen-containing functional groups
            if c_count == 1 && h_count == 4 && elements["O"] == 1
                # Methanol (CH3OH)
                # Replace default values with more accurate group contributions
                cp_298 = group_contrib["C-(C)(H)3"][1] + group_contrib["O-(C)(H)"][1]
                h_298 = group_contrib["C-(C)(H)3"][2] + group_contrib["O-(C)(H)"][2]
                s_298 = group_contrib["C-(C)(H)3"][3] + group_contrib["O-(C)(H)"][3]
                
                # Update Cp coefficients
                a_sum = cp_coeffs["C-(C)(H)3"][1] + cp_coeffs["O-(C)(H)"][1]
                b_sum = cp_coeffs["C-(C)(H)3"][2] + cp_coeffs["O-(C)(H)"][2]
                c_sum = cp_coeffs["C-(C)(H)3"][3] + cp_coeffs["O-(C)(H)"][3]
                d_sum = cp_coeffs["C-(C)(H)3"][4] + cp_coeffs["O-(C)(H)"][4]
            end
        end
    end
    
    # Calculate Cp at the requested temperature using the polynomial
    cp = a_sum + b_sum * temperature + c_sum * temperature^2 + d_sum * temperature^3
    
    # If the calculated Cp is unrealistic, use the 298K value with a temperature correction
    if cp <= 0 || cp > 1000
        cp = cp_298 * (1.0 + 0.0015 * (temperature - 298.15))
    end
    
    # Calculate enthalpy change from 298.15 K to T
    h_change = a_sum * (temperature - 298.15) +
               b_sum * (temperature^2 - 298.15^2) / 2 +
               c_sum * (temperature^3 - 298.15^3) / 3 +
               d_sum * (temperature^4 - 298.15^4) / 4
    
    # Calculate total enthalpy
    h = h_298 + h_change
    
    # Calculate entropy change from 298.15 K to T
    s_change = a_sum * log(temperature / 298.15) +
               b_sum * (temperature - 298.15) +
               c_sum * (temperature^2 - 298.15^2) / 2 +
               d_sum * (temperature^3 - 298.15^3) / 3
    
    # Calculate total entropy
    s = s_298 + s_change
    
    # Calculate Gibbs free energy
    g = h - temperature * s / 1000  # Convert S from J/mol/K to kJ/mol/K
    
    # Estimate uncertainties - Benson method has better accuracy than simple methods
    cp_uncertainty = max(cp * 0.05, 2.0)  # 5% or at least 2 J/mol/K
    h_uncertainty = max(abs(h) * 0.08, 5.0)  # 8% or at least 5 kJ/mol
    s_uncertainty = max(s * 0.04, 4.0)  # 4% or at least 4 J/mol/K
    g_uncertainty = max(abs(g) * 0.08, 5.0)  # 8% or at least 5 kJ/mol
    
    # Create result with uncertainties
    return Dict(
        "temperature" => temperature,
        "Cp" => cp,
        "Cp_uncertainty" => cp_uncertainty,
        "H" => h,
        "H_uncertainty" => h_uncertainty,
        "S" => s,
        "S_uncertainty" => s_uncertainty,
        "G" => g,
        "G_uncertainty" => g_uncertainty,
        "method" => "benson_group_additivity"
    )
end

"""
    joback_method(formula::String, temperature::Float64)

Estimate thermodynamic properties using the Joback group contribution method.
This is a more accurate method for organic compounds when formula can be parsed into groups.
"""
function joback_method(formula::String, temperature::Float64)
    # This function attempts to parse the formula into Joback groups
    # and calculate thermodynamic properties
    
    # Parse formula to identify elements
    elements = extract_formula_elements(formula)
    
    if length(elements) == 0
        error("Could not parse formula: $formula")
    end
    
    # Define Joback group contributions
    # These are the original Joback parameter values for various functional groups
    joback_groups = Dict(
        # Group name => [ΔHf, ΔGf, a, b, c, d]
        # where Cp = a + bT + cT^2 + dT^3
        # Units: ΔHf, ΔGf in kJ/mol; a, b, c, d for Cp in J/mol/K
        
        # Non-ring groups
        "-CH3" => [      -76.45,    -43.96,    19.5,    -8.08e-3,   1.53e-4,   -9.67e-8],
        "-CH2-" => [     -20.64,     8.42,     -0.909,   9.50e-2,  -5.44e-5,    1.19e-8],
        ">CH-" => [       29.89,    58.36,     -23.0,    2.04e-1,  -2.65e-4,    1.20e-7],
        ">C<" => [        82.23,   116.02,     -66.2,    4.27e-1,  -6.41e-4,    3.01e-7],
        "=CH2" => [       -9.63,     3.77,     23.6,     -3.81e-2,   1.72e-4,   -1.03e-7],
        "=CH-" => [       37.97,    48.53,     -8.00,    1.05e-1,  -9.63e-5,    3.56e-8],
        "=C<" => [        83.99,    92.36,     -8.83,    2.48e-2,   4.46e-5,   -6.32e-8],
        "≡CH" => [       132.13,   125.45,     24.5,    -2.71e-2,   1.11e-4,   -6.78e-8],
        "≡C-" => [       142.14,   146.70,     24.5,    -2.71e-2,   1.11e-4,   -6.78e-8],
        
        # Oxygen groups
        "-OH (alcohol)" => [-208.04, -189.20,   25.7,    -6.91e-2,   1.77e-4,   -9.88e-8],
        "-OH (phenol)" => [-221.65, -197.37,   -2.81,    1.11e-1,  -1.16e-4,    4.94e-8],
        "-O- (non-ring)" => [-132.22, -105.00,  22.2,    -9.09e-2,   1.22e-4,   -4.29e-8],
        "-O- (ring)" => [  -138.16, -110.50,   10.5,    -1.05e-1,   2.03e-4,   -1.31e-7],
        ">C=O (non-ring)" => [-133.22, -120.50, 6.45,     6.70e-2,  -3.57e-5,    2.86e-9],
        ">C=O (ring)" => [ -164.50, -126.27,   31.1,    -1.32e-2,   1.32e-4,   -9.86e-8],
        "O=CH- (aldehyde)" => [-162.03, -143.48, 30.9, -3.36e-2,  1.60e-4,   -1.03e-7],
        "-COOH (acid)" => [-426.72, -387.87,   24.1,     4.27e-2,   8.04e-5,   -6.87e-8],
        "-COO- (ester)" => [-337.92, -301.95,   24.5,     4.02e-2,   4.02e-5,   -4.52e-8],
        
        # Nitrogen groups
        "-NH2" => [       -22.02,    38.81,     26.9,    -4.12e-2,   1.64e-4,   -9.76e-8],
        ">NH (non-ring)" => [73.23,   93.70,     5.69,     7.31e-2,  -4.72e-5,    1.45e-8],
        ">NH (ring)" => [   50.17,    52.66,     7.24,     2.45e-2,   3.70e-5,   -4.12e-8],
        ">N- (non-ring)" => [123.34, 163.16,   -2.03,     1.32e-1,  -1.41e-4,    6.20e-8],
        "-N= (non-ring)" => [23.61,    23.91,     33.5,    -2.99e-2,   1.17e-4,   -7.52e-8],
        "-CN" => [        125.66,   139.12,     33.6,     2.35e-2,   3.80e-5,   -4.15e-8],
        "-NO2" => [       -66.57,   -16.83,     25.5,     4.13e-2,   3.05e-5,   -3.67e-8],
        
        # Sulfur groups
        "-SH" => [        -17.33,    -8.48,     28.5,    -4.16e-2,   1.19e-4,   -6.37e-8],
        "-S- (non-ring)" => [41.87,   33.12,     16.6,     1.95e-2,   1.89e-5,   -2.85e-8],
        
        # Halogen groups
        "-F" => [        -251.92,  -247.19,     27.1,    -1.34e-2,   1.18e-4,   -8.39e-8],
        "-Cl" => [       -71.55,   -64.31,     36.0,    -7.33e-2,   1.84e-4,   -1.03e-7],
        "-Br" => [       -29.48,   -38.06,     42.9,    -7.63e-2,   1.82e-4,   -1.04e-7],
        "-I" => [         21.06,    14.15,     41.7,    -6.87e-2,   1.67e-4,   -9.88e-8]
    )
    
    # Try to identify groups based on formula
    # This is a very simplified approach - a real implementation would use
    # a proper formula parser or SMILES representation
    
    # Initialize group counts
    group_counts = Dict{String, Int}()
    
    # For simple molecules, try to identify known structures
    if haskey(elements, "C") && haskey(elements, "H")
        c_count = elements["C"]
        h_count = elements["H"]
        
        if c_count == 1 && h_count == 4
            # Methane (CH4)
            group_counts["-CH3"] = 1
        elseif c_count == 2 && h_count == 6
            # Ethane (C2H6)
            group_counts["-CH3"] = 2
        elseif c_count == 2 && h_count == 4
            # Ethylene (C2H4)
            group_counts["=CH2"] = 2
        elseif c_count == 2 && h_count == 2
            # Acetylene (C2H2)
            group_counts["≡CH"] = 2
        elseif c_count == 3 && h_count == 8
            # Propane (C3H8)
            group_counts["-CH3"] = 2
            group_counts["-CH2-"] = 1
        elseif c_count >= 4 && h_count == 2 * c_count + 2
            # Longer alkanes: CnH(2n+2)
            group_counts["-CH3"] = 2
            group_counts["-CH2-"] = c_count - 2
        elseif haskey(elements, "O") && c_count == 1 && h_count == 4 && elements["O"] == 1
            # Methanol (CH3OH)
            group_counts["-CH3"] = 1
            group_counts["-OH (alcohol)"] = 1
        elseif haskey(elements, "O") && c_count == 2 && h_count == 6 && elements["O"] == 1
            # Ethanol (C2H5OH)
            group_counts["-CH3"] = 1
            group_counts["-CH2-"] = 1
            group_counts["-OH (alcohol)"] = 1
        elseif haskey(elements, "O") && c_count == 1 && h_count == 2 && elements["O"] == 1
            # Formaldehyde (CH2O)
            group_counts["O=CH- (aldehyde)"] = 1
        elseif haskey(elements, "O") && c_count == 2 && h_count == 4 && elements["O"] == 1
            # Acetaldehyde (CH3CHO)
            group_counts["-CH3"] = 1
            group_counts["O=CH- (aldehyde)"] = 1
        elseif haskey(elements, "O") && c_count == 2 && h_count == 4 && elements["O"] == 2
            # Acetic acid (CH3COOH)
            group_counts["-CH3"] = 1
            group_counts["-COOH (acid)"] = 1
        else
            # For more complex molecules, use element-based estimation
            # Estimate group counts based on element counts
            if c_count > 0
                # Assume most carbons are in methylene groups
                group_counts["-CH2-"] = max(0, c_count - 2)
                
                # Add two methyl end groups if enough hydrogens
                if h_count >= 6 && c_count >= 2
                    group_counts["-CH3"] = 2
                    group_counts["-CH2-"] = max(0, c_count - 2)
                else
                    # Not enough H for all CH2, so add some CH groups
                    group_counts["-CH2-"] = max(0, min(group_counts["-CH2-"], h_count ÷ 2))
                    group_counts[">CH-"] = max(0, min(c_count - group_counts["-CH2-"], h_count - 2 * group_counts["-CH2-"]))
                    group_counts[">C<"] = max(0, c_count - group_counts["-CH2-"] - group_counts[">CH-"])
                end
            end
            
            # Add functional groups for other elements
            if haskey(elements, "O")
                o_count = elements["O"]
                if o_count == 1
                    # One oxygen - could be alcohol, ether, or carbonyl
                    if h_count >= c_count * 2  # Enough H for alcohol
                        group_counts["-OH (alcohol)"] = 1
                    else  # Probably carbonyl
                        group_counts[">C=O (non-ring)"] = 1
                    end
                elseif o_count == 2
                    # Two oxygens - could be acid, ester, or two other groups
                    if c_count >= 2 && h_count <= c_count * 2  # Likely acid
                        group_counts["-COOH (acid)"] = 1
                    else  # Maybe two alcohols
                        group_counts["-OH (alcohol)"] = min(2, o_count)
                    end
                else
                    # Multiple oxygens - add a mix of groups
                    group_counts["-OH (alcohol)"] = min(2, o_count)
                    if o_count > 2
                        group_counts["-O- (non-ring)"] = min(1, o_count - 2)
                        if o_count > 3
                            group_counts[">C=O (non-ring)"] = min(1, o_count - 3)
                        end
                    end
                end
            end
            
            # Add nitrogen groups
            if haskey(elements, "N")
                n_count = elements["N"]
                if n_count == 1
                    # One nitrogen - could be amine, nitrile, etc.
                    if c_count >= 1 && h_count >= 2  # Likely amine
                        group_counts["-NH2"] = 1
                    else if c_count >= 1  # Maybe nitrile
                        group_counts["-CN"] = 1
                    end
                else
                    # Multiple nitrogens - add a mix of groups
                    group_counts["-NH2"] = min(1, n_count)
                    if n_count > 1
                        group_counts[">NH (non-ring)"] = min(1, n_count - 1)
                    end
                end
            end
            
            # Add halogen groups
            if haskey(elements, "F")
                group_counts["-F"] = elements["F"]
            end
            if haskey(elements, "Cl")
                group_counts["-Cl"] = elements["Cl"]
            end
            if haskey(elements, "Br")
                group_counts["-Br"] = elements["Br"]
            end
            if haskey(elements, "I")
                group_counts["-I"] = elements["I"]
            end
        end
    else
        # Non-hydrocarbon molecule
        # Use element-based approach for less common molecules
        if haskey(elements, "C")
            group_counts["-CH3"] = min(2, elements["C"])
            if elements["C"] > 2
                group_counts["-CH2-"] = elements["C"] - 2
            end
        end
        
        if haskey(elements, "O")
            group_counts["-OH (alcohol)"] = min(2, elements["O"])
        end
        
        if haskey(elements, "N")
            group_counts["-NH2"] = min(1, elements["N"])
        end
        
        if haskey(elements, "S")
            group_counts["-SH"] = min(1, elements["S"])
        end
        
        if haskey(elements, "F")
            group_counts["-F"] = elements["F"]
        end
        
        if haskey(elements, "Cl")
            group_counts["-Cl"] = elements["Cl"]
        end
        
        if haskey(elements, "Br")
            group_counts["-Br"] = elements["Br"]
        end
        
        if haskey(elements, "I")
            group_counts["-I"] = elements["I"]
        end
    end
    
    # If no groups were identified, can't use Joback method
    if isempty(group_counts)
        error("Could not identify Joback groups in formula: $formula")
    end
    
    # Sum the contributions from each group
    delta_hf_sum = 0.0
    delta_gf_sum = 0.0
    a_sum = 0.0
    b_sum = 0.0
    c_sum = 0.0
    d_sum = 0.0
    
    for (group, count) in group_counts
        if count > 0 && haskey(joback_groups, group)
            params = joback_groups[group]
            delta_hf_sum += params[1] * count
            delta_gf_sum += params[2] * count
            a_sum += params[3] * count
            b_sum += params[4] * count
            c_sum += params[5] * count
            d_sum += params[6] * count
        end
    end
    
    # Joback equations for estimating properties
    # Standard heat of formation (kJ/mol)
    delta_hf = 68.29 + delta_hf_sum
    
    # Standard Gibbs energy of formation (kJ/mol)
    delta_gf = 53.88 + delta_gf_sum
    
    # Calculate heat capacity at the given temperature
    # Cp = a + bT + cT^2 + dT^3
    t_kelvin = temperature  # K
    cp = (a_sum + b_sum * t_kelvin + c_sum * t_kelvin^2 + d_sum * t_kelvin^3)
    
    # Standard entropy at 298.15 K
    # S298 = 198.2 + ∑(contributions)
    # For simplicity, using a correlation based on functional groups
    s_298 = 198.2 + a_sum / 2  # Simplified correlation
    
    # Calculate entropy at the given temperature
    # For accurate calculation, we would need to integrate Cp/T from 298.15 to T
    if t_kelvin != 298.15
        # Simplified temperature correction
        s_t = s_298 + cp * log(t_kelvin / 298.15)
    else
        s_t = s_298
    end
    
    # Calculate enthalpy at temperature T
    h_298 = delta_hf
    h_t = h_298 + a_sum * (t_kelvin - 298.15) + 
                  b_sum * (t_kelvin^2 - 298.15^2) / 2 + 
                  c_sum * (t_kelvin^3 - 298.15^3) / 3 + 
                  d_sum * (t_kelvin^4 - 298.15^4) / 4
    
    # Calculate Gibbs free energy at temperature T
    g_t = h_t - t_kelvin * s_t / 1000  # Convert S from J/mol/K to kJ/mol/K
    
    # Estimate uncertainties based on the number of functional groups
    # More groups generally leads to higher uncertainty
    num_groups = sum(values(group_counts))
    cp_uncertainty = max(cp * 0.08, 3.0) * (1 + 0.03 * num_groups)
    h_uncertainty = max(abs(h_t) * 0.1, 8.0) * (1 + 0.03 * num_groups)
    s_uncertainty = max(s_t * 0.08, 5.0) * (1 + 0.02 * num_groups)
    g_uncertainty = max(abs(g_t) * 0.1, 8.0) * (1 + 0.03 * num_groups)
    
    # Create result with uncertainties
    return Dict(
        "temperature" => temperature,
        "Cp" => cp,
        "Cp_uncertainty" => cp_uncertainty,
        "H" => h_t,
        "H_uncertainty" => h_uncertainty,
        "S" => s_t,
        "S_uncertainty" => s_uncertainty,
        "G" => g_t,
        "G_uncertainty" => g_uncertainty,
        "method" => "joback_method",
        "identified_groups" => group_counts
    )
end