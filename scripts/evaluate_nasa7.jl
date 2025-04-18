#!/usr/bin/env julia

# This script demonstrates how to evaluate NASA-7 polynomial for N2 
# without relying on complex database queries

using Pkg
Pkg.activate(".")

println("NASA-7 Polynomial Evaluation for N2")
println("===================================")

# Physical constants
const R = 8.31446261815324  # Universal gas constant, J/(mol·K)

# NASA-7 polynomial coefficients for N2 from GRI-Mech 3.0
# Range 1: 200-1000K
# Range 2: 1000-6000K
nasa_coeff = [
    # Low temperature range (200-1000K)
    [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372],
    # High temperature range (1000-6000K)
    [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]
]

temperature_ranges = [[200.0, 1000.0], [1000.0, 6000.0]]

# Function to calculate Cp/R
function calculate_nasa7_cp(coeffs, T)
    return coeffs[1] + coeffs[2]*T + coeffs[3]*T^2 + coeffs[4]*T^3 + coeffs[5]*T^4
end

# Function to calculate H/RT
function calculate_nasa7_enthalpy(coeffs, T)
    return coeffs[1] + (coeffs[2]/2)*T + (coeffs[3]/3)*T^2 + (coeffs[4]/4)*T^3 + (coeffs[5]/5)*T^4 + coeffs[6]/T
end

# Function to calculate S/R
function calculate_nasa7_entropy(coeffs, T)
    return coeffs[1]*log(T) + coeffs[2]*T + (coeffs[3]/2)*T^2 + (coeffs[4]/3)*T^3 + (coeffs[5]/4)*T^4 + coeffs[7]
end

# Function to get coefficients for a specific temperature
function get_coeffs_for_temperature(coeff_sets, temp_ranges, T)
    for i in 1:length(temp_ranges)
        if temp_ranges[i][1] <= T && T <= temp_ranges[i][2]
            return coeff_sets[i]
        end
    end
    error("Temperature $T is outside all defined ranges")
end

# Temperature to evaluate
T = 298.15  # K

# Get appropriate coefficients
coeffs = get_coeffs_for_temperature(nasa_coeff, temperature_ranges, T)

# Calculate thermodynamic properties
cp_over_R = calculate_nasa7_cp(coeffs, T)
h_over_RT = calculate_nasa7_enthalpy(coeffs, T)
s_over_R = calculate_nasa7_entropy(coeffs, T)

# Convert to standard units
cp = cp_over_R * R             # J/mol/K
h = h_over_RT * R * T / 1000   # kJ/mol
s = s_over_R * R               # J/mol/K
g = h - T * s / 1000           # kJ/mol

# Display results
println("Properties for N2 at $T K (NASA-7 polynomial):")
println("Cp = $cp J/mol/K")
println("H = $h kJ/mol")
println("S = $s J/mol/K")
println("G = $g kJ/mol")

# Compare to reference values (from NIST)
println("\nComparison with reference values (NIST):")
println("Cp (NIST) = 29.12 J/mol/K")
println("S (NIST) = 191.61 J/mol/K")
println("H (NIST) = 8.67 kJ/mol (relative to 0K)")
println("\nNote: Our enthalpy calculations may differ from NIST as they use a different reference state.")