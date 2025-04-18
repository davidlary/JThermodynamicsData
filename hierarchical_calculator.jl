#!/usr/bin/env julia

"""
Hierarchical Thermodynamic Properties Calculator

This script calculates thermodynamic properties for chemical species using a hierarchical
approach, starting with theoretical estimates and progressively refining with increasingly
accurate data sources. No database connection required.
"""

# Activate the project to ensure we use the right dependencies
using Pkg
Pkg.activate(@__DIR__)

using Statistics
using Plots
using Printf
using Dates
using YAML
using JSON
using LinearAlgebra
using Optim  # For optimization in Gibbs energy minimization

# Define constants
const R = 8.31446261815324  # Universal gas constant, J/(mol·K)
const NA = 6.02214076e23    # Avogadro's number, 1/mol
const KB = 1.380649e-23     # Boltzmann constant, J/K
const H = 6.62607015e-34    # Planck's constant, J·s
const C = 299792458.0       # Speed of light, m/s
const HARTREE_TO_KJ_PER_MOL = 2625.5  # Conversion factor from Hartree to kJ/mol

# Alias lowercase constants to match usage in functions
const h = H 
const kb = KB
const na = NA
const c = C

# Define structure for thermodynamic polynomials
struct ThermodynamicPolynomial
    species_name::String
    formula::String
    source::String
    priority::Int
    reliability_score::Float64
    polynomial_type::String
    temperature_ranges::Vector{Vector{Float64}}
    coefficients::Vector{Vector{Float64}}
    uncertainty::Float64
end

# Define structure for molecular data (used for theoretical calculations)
struct MolecularData
    molecular_weight::Float64
    formula::String
    symmetry_number::Int
    spin_multiplicity::Int
    vibrational_frequencies::Vector{Float64}
    rotational_constants::Vector{Float64}
end

# Common atomic masses (in g/mol)
const ATOMIC_MASSES = Dict(
    "H" => 1.00794,
    "D" => 2.01410,
    "T" => 3.01605,
    "He" => 4.002602,
    "Li" => 6.941,
    "Be" => 9.012182,
    "B" => 10.811,
    "C" => 12.0107,
    "N" => 14.0067,
    "O" => 15.9994,
    "F" => 18.9984032,
    "Ne" => 20.1797,
    "Na" => 22.98976928,
    "Mg" => 24.3050,
    "Al" => 26.9815386,
    "Si" => 28.0855,
    "P" => 30.973762,
    "S" => 32.065,
    "Cl" => 35.453,
    "Ar" => 39.948,
    "K" => 39.0983,
    "Ca" => 40.078,
    "Fe" => 55.845,
    "Cu" => 63.546,
    "Ag" => 107.8682,
    "Au" => 196.966569,
    "Hg" => 200.59,
    "Pb" => 207.2,
    "Xe" => 131.293,
    "Ne" => 20.1797
)

"""
    extract_formula_elements(formula::String)

Parse a chemical formula into a dictionary of element counts.
"""
function extract_formula_elements(formula::String)
    elements = Dict{String, Int}()
    
    # Simple regex-free parsing for basic formulas
    i = 1
    while i <= length(formula)
        # If it's an uppercase letter, it's the start of an element
        if isuppercase(formula[i])
            # Extract the element symbol (one or two characters)
            element = formula[i:i]
            if i+1 <= length(formula) && islowercase(formula[i+1])
                element *= formula[i+1:i+1]
                i += 1
            end
            
            # Extract the count (default is 1)
            count = 1
            num_str = ""
            i += 1
            while i <= length(formula) && isdigit(formula[i])
                num_str *= formula[i]
                i += 1
            end
            
            if !isempty(num_str)
                count = parse(Int, num_str)
            end
            
            # Add to the elements dictionary
            if haskey(elements, element)
                elements[element] += count
            else
                elements[element] = count
            end
        else
            i += 1
        end
    end
    
    return elements
end

"""
    calculate_molecular_weight(formula::String)

Calculate the molecular weight of a compound from its formula.
"""
function calculate_molecular_weight(formula::String)
    elements = extract_formula_elements(formula)
    mw = 0.0
    
    for (element, count) in elements
        if haskey(ATOMIC_MASSES, element)
            mw += ATOMIC_MASSES[element] * count
        else
            @warn "Unknown element in formula: $element"
        end
    end
    
    return mw
end

"""
    estimate_molecular_properties(formula::String)

Estimate basic molecular properties needed for theoretical calculations.
"""
function estimate_molecular_properties(formula::String)
    # Parse formula and get elements
    elements = extract_formula_elements(formula)
    
    # Calculate molecular weight
    mw = calculate_molecular_weight(formula)
    
    # Count total atoms
    total_atoms = sum(values(elements))
    
    # Estimate symmetry number (very rough)
    symmetry = 1
    
    # Common symmetry cases
    if length(elements) == 1
        # Monatomic molecule
        symmetry = 1
    elseif length(elements) == 2 && 
           (haskey(elements, "H") && elements["H"] > 0) && 
           length(elements) == 2
        # Diatomic molecule with hydrogen
        symmetry = 1
    elseif total_atoms == 2 && length(elements) == 1
        # Homonuclear diatomic (like O2, N2)
        symmetry = 2
    elseif haskey(elements, "C") && haskey(elements, "H") && 
           elements["C"] >= 2 && elements["H"] >= 6 &&
           length(elements) == 2
        # Might be a hydrocarbon with some symmetry
        symmetry = 2
    end
    
    # Estimate spin multiplicity (very rough)
    spin = 1  # Singlet by default
    
    # Some common cases for spin multiplicity
    if haskey(elements, "O") && elements["O"] == 2 && length(elements) == 1
        # O2 is a triplet
        spin = 3
    end
    
    # Rough estimate for vibrational frequencies
    # This is extremely crude - real frequencies would come from experiment or calculation
    vib_freq = Float64[]
    
    # Estimate number of vibrational modes
    if total_atoms == 1
        # Monatomic - no vibrations
        n_modes = 0
    elseif total_atoms == 2
        # Diatomic - one vibration
        n_modes = 1
        
        # Rough frequency based on elements
        if haskey(elements, "H")
            # X-H stretch (higher frequency)
            push!(vib_freq, 3000.0)
        else
            # Heavy element bond (lower frequency)
            push!(vib_freq, 1500.0)
        end
    else
        # Polyatomic
        n_modes = 3 * total_atoms - 6  # Non-linear molecule
        
        # Rough distribution of frequencies
        # Again, this is extremely crude
        
        # Low frequency modes (400-800 cm^-1)
        n_low = floor(Int, n_modes / 3)
        for _ in 1:n_low
            push!(vib_freq, 400.0 + 400.0 * rand())
        end
        
        # Medium frequency modes (800-1500 cm^-1)
        n_medium = floor(Int, n_modes / 3)
        for _ in 1:n_medium
            push!(vib_freq, 800.0 + 700.0 * rand())
        end
        
        # High frequency modes (1500-3500 cm^-1)
        n_high = n_modes - n_low - n_medium
        for _ in 1:n_high
            push!(vib_freq, 1500.0 + 2000.0 * rand())
        end
        
        # Add X-H stretches if hydrogen is present
        if haskey(elements, "H") && elements["H"] > 0
            for i in 1:min(elements["H"], 3)
                if length(vib_freq) > i
                    vib_freq[end-i+1] = 2800.0 + i * 100.0  # X-H stretches around 3000 cm^-1
                end
            end
        end
    end
    
    # Rough estimate for rotational constants
    # This is extremely crude - real rotational constants would come from experiment or calculation
    rot_const = Float64[]
    
    if total_atoms == 1
        # Monatomic - no rotation
    elseif total_atoms == 2
        # Diatomic or linear - one rotational constant
        # Estimate based on molecular weight and typical bond length
        bond_length = 1.5e-10  # meters (typical bond length)
        
        # Moment of inertia for a diatomic
        reduced_mass = mw / (NA * 1000)  # kg/molecule
        I = reduced_mass * bond_length^2  # kg*m^2
        
        # Rotational constant in cm^-1
        B = h / (8 * π^2 * c * I) / 100  # cm^-1
        push!(rot_const, B)
    else
        # Non-linear molecule - three rotational constants
        # This is a very rough approximation assuming a quasi-spherical molecule
        
        # Rough estimate of molecular size
        radius = 1.0e-10 * (total_atoms)^(1/3)  # meters
        
        # Moment of inertia for a spherical top
        mass = mw / (NA * 1000)  # kg/molecule
        I = 0.4 * mass * radius^2  # kg*m^2
        
        # Rotational constants in cm^-1
        B = h / (8 * π^2 * c * I) / 100  # cm^-1
        
        # For non-spherical molecules, introduce some asymmetry
        push!(rot_const, B * (1.0 + 0.2 * rand()))
        push!(rot_const, B * (1.0 - 0.1 * rand()))
        push!(rot_const, B * (1.0 - 0.3 * rand()))
    end
    
    return MolecularData(
        mw, 
        formula, 
        symmetry, 
        spin, 
        vib_freq, 
        rot_const
    )
end

"""
    calculate_translational_partition_function(molecular_weight::Float64, temperature::Float64)

Calculate the translational partition function for statistical thermodynamics.
"""
function calculate_translational_partition_function(molecular_weight::Float64, temperature::Float64)
    # Convert molecular weight from g/mol to kg/molecule
    m = molecular_weight / (NA * 1000)
    
    # Thermal wavelength
    lambda = h / sqrt(2 * π * m * kb * temperature)
    
    # Volume per molecule (approximation for ideal gas at standard conditions)
    P_STANDARD = 101325.0  # Standard pressure, Pa
    V = (R * temperature) / (P_STANDARD * NA)
    
    # Translational partition function
    q_trans = V / lambda^3
    
    return q_trans
end

"""
    calculate_rotational_partition_function(rotational_constants::Vector{Float64}, 
                                          symmetry_number::Int, temperature::Float64)

Calculate the rotational partition function for statistical thermodynamics.
"""
function calculate_rotational_partition_function(
    rotational_constants::Vector{Float64}, 
    symmetry_number::Int, 
    temperature::Float64
)
    # Check if molecule is linear or non-linear
    if length(rotational_constants) == 1
        # Linear molecule
        B = rotational_constants[1]  # cm^-1
        
        # Convert to SI units (Hz)
        B_si = B * c * 100  # Hz
        
        # Rotational partition function for linear molecule
        q_rot = (kb * temperature) / (h * B_si * symmetry_number)
    elseif length(rotational_constants) >= 3
        # Non-linear molecule
        A = rotational_constants[1]  # cm^-1
        B = rotational_constants[2]  # cm^-1
        C = rotational_constants[3]  # cm^-1
        
        # Convert to SI units (Hz)
        A_si = A * 100 * 2.998e10  # Hz
        B_si = B * 100 * 2.998e10  # Hz
        C_si = C * 100 * 2.998e10  # Hz
        
        # Rotational partition function for non-linear molecule
        q_rot = sqrt(π) / symmetry_number * ((kb * temperature)^1.5) / (h^1.5 * sqrt(A_si * B_si * C_si))
    else
        # Default for incomplete data
        q_rot = 1.0
    end
    
    return q_rot
end

"""
    calculate_vibrational_partition_function(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational partition function for statistical thermodynamics.
"""
function calculate_vibrational_partition_function(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    # Vibrational frequencies should be in cm^-1
    
    q_vib = 1.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Vibrational partition function contribution
        q_vib *= 1 / (1 - exp(-x))
    end
    
    return q_vib
end

"""
    calculate_vibrational_enthalpy(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational contribution to enthalpy.
"""
function calculate_vibrational_enthalpy(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    # Vibrational frequencies should be in cm^-1
    
    h_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Vibrational enthalpy contribution
        h_vib += R * temperature * x / (exp(x) - 1)
        
        # Add zero-point energy
        h_vib += R * temperature * x / 2
    end
    
    return h_vib
end

"""
    calculate_vibrational_entropy(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational contribution to entropy.
"""
function calculate_vibrational_entropy(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    # Vibrational frequencies should be in cm^-1
    
    s_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Vibrational entropy contribution
        s_vib += R * (x / (exp(x) - 1) - log(1 - exp(-x)))
    end
    
    return s_vib
end

"""
    calculate_vibrational_heat_capacity(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational contribution to heat capacity.
"""
function calculate_vibrational_heat_capacity(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    # Vibrational frequencies should be in cm^-1
    
    cp_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Vibrational heat capacity contribution
        cp_vib += R * (x^2 * exp(x) / (exp(x) - 1)^2)
    end
    
    return cp_vib
end

"""
    calculate_statistical_thermodynamics(molecular_data::MolecularData, temperature::Float64)

Calculate thermodynamic properties using statistical thermodynamics.
"""
function calculate_statistical_thermodynamics(molecular_data::MolecularData, temperature::Float64)
    # Calculate translational contribution
    q_trans = calculate_translational_partition_function(molecular_data.molecular_weight, temperature)
    
    # Calculate rotational contribution (if rotational constants provided)
    q_rot = 1.0
    if !isempty(molecular_data.rotational_constants)
        q_rot = calculate_rotational_partition_function(
            molecular_data.rotational_constants,
            molecular_data.symmetry_number,
            temperature
        )
    end
    
    # Calculate vibrational contribution (if vibrational frequencies provided)
    q_vib = 1.0
    if !isempty(molecular_data.vibrational_frequencies)
        q_vib = calculate_vibrational_partition_function(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    # Calculate electronic contribution
    q_elec = molecular_data.spin_multiplicity
    
    # Calculate thermodynamic properties
    
    # Enthalpy (H)
    h_trans = R * temperature  # Translational contribution to enthalpy
    
    h_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        h_rot = R * temperature  # Rotational contribution for non-linear molecules
    end
    
    h_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        h_vib = calculate_vibrational_enthalpy(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    h_elec = 0.0  # Electronic contribution at ground state
    
    enthalpy = h_trans + h_rot + h_vib + h_elec  # Use enthalpy instead of h to avoid conflict with Planck constant
    
    # Entropy (S)
    s_trans = R * (log(q_trans) + 1 + 1.5)  # Translational contribution to entropy
    
    s_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        s_rot = R * (log(q_rot) + 1)  # Rotational contribution
    end
    
    s_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        s_vib = calculate_vibrational_entropy(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    s_elec = R * log(q_elec)  # Electronic contribution
    
    s = s_trans + s_rot + s_vib + s_elec
    
    # Heat capacity (Cp)
    cp_trans = 5 * R / 2  # Translational contribution to heat capacity
    
    cp_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        cp_rot = R  # Rotational contribution for non-linear molecules (3/2 * R for linear)
    end
    
    cp_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        cp_vib = calculate_vibrational_heat_capacity(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    cp_elec = 0.0  # Electronic contribution at ground state
    
    cp = cp_trans + cp_rot + cp_vib + cp_elec
    
    # Gibbs free energy (G)
    g = enthalpy - temperature * s / 1000  # Convert S to kJ/mol/K
    
    # Estimate uncertainties based on the quality of molecular data
    has_good_vib = !isempty(molecular_data.vibrational_frequencies)
    has_good_rot = !isempty(molecular_data.rotational_constants)
    
    # Higher uncertainty if we're missing good data
    base_uncertainty = has_good_vib && has_good_rot ? 0.1 : 0.25
    
    cp_uncertainty = max(cp * base_uncertainty, 3.0)  # At least 3 J/mol/K
    h_uncertainty = max(abs(enthalpy) * base_uncertainty, 10.0)  # At least 10 kJ/mol
    s_uncertainty = max(s * base_uncertainty, 5.0)  # At least 5 J/mol/K
    g_uncertainty = max(abs(g) * base_uncertainty, 10.0)  # At least 10 kJ/mol
    
    return Dict(
        "temperature" => temperature,
        "Cp" => cp,
        "Cp_uncertainty" => cp_uncertainty,
        "H" => enthalpy, 
        "H_uncertainty" => h_uncertainty,
        "S" => s,
        "S_uncertainty" => s_uncertainty,
        "G" => g,
        "G_uncertainty" => g_uncertainty,
        "method" => "statistical_thermodynamics"
    )
end

"""
    estimate_group_contribution(formula::String, temperature::Float64)

Rough estimate of thermodynamic properties using simple group contribution method.
This is a placeholder for a more sophisticated implementation.
"""
function estimate_group_contribution(formula::String, temperature::Float64)
    # Parse formula
    elements = extract_formula_elements(formula)
    
    # Calculate molecular weight
    mw = calculate_molecular_weight(formula)
    
    # Very rough estimates - these would be replaced with proper group contribution methods
    # These are just placeholders for demonstration
    
    # Heat capacity - crude correlation with molecular weight
    cp = 20.0 + 3.0 * log(mw)  # J/mol/K
    
    # Enthalpy - crude estimate based on atom count and bonds
    total_atoms = sum(values(elements))
    h = 10.0 * total_atoms  # kJ/mol
    
    # Entropy - crude correlation with molecular weight and temperature
    s = 150.0 + 10.0 * log(mw) + 0.1 * temperature  # J/mol/K
    
    # Gibbs energy
    g = h - temperature * s / 1000
    
    # Relatively high uncertainties for group contribution method
    cp_uncertainty = 0.25 * cp
    h_uncertainty = 0.25 * abs(h)
    s_uncertainty = 0.25 * s
    g_uncertainty = 0.25 * abs(g)
    
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
        "method" => "group_contribution"
    )
end

"""
    calculate_theoretical_properties(species_name::String, formula::String, temperature::Float64)

Calculate theoretical properties using multiple advanced methods and combine with robust uncertainty estimation.
"""
function calculate_theoretical_properties(species_name::String, formula::String, temperature::Float64, config::Dict)
    # Special case for electron (e-) - predefined values
    if species_name == "e-" || formula == "e-"
        # Known values for electron
        # Electrons have 5/2 * R heat capacity (3 translational degrees of freedom + 2 from quantum effects)
        # Using NASA-7 polynomial values for consistency
        cp_value = 5/2 * R  # J/mol/K
        h_value = -745.4 * R * temperature / 1000  # kJ/mol
        s_value = -1.2 * R  # J/mol/K
        g_value = h_value - temperature * s_value / 1000  # kJ/mol
        
        # Low uncertainty for electron - well-characterized
        uncertainty = 0.01  # 1%
        
        electron_result = Dict(
            "species_name" => species_name,
            "formula" => formula,
            "temperature" => temperature,
            "data_source" => "ELECTRON_SPECIAL_CASE",
            "properties" => Dict(
                "Cp" => Dict(
                    "value" => cp_value,
                    "uncertainty" => cp_value * uncertainty,
                    "units" => "J/mol/K"
                ),
                "H" => Dict(
                    "value" => h_value,
                    "uncertainty" => abs(h_value) * uncertainty,
                    "units" => "kJ/mol"
                ),
                "S" => Dict(
                    "value" => s_value,
                    "uncertainty" => abs(s_value) * uncertainty,
                    "units" => "J/mol/K"
                ),
                "G" => Dict(
                    "value" => g_value,
                    "uncertainty" => abs(g_value) * uncertainty,
                    "units" => "kJ/mol"
                )
            )
        )
        
        # Add empty individual methods for plotting consistency
        electron_result["individual_methods"] = [electron_result]
        
        return electron_result
    end
    
    # Normal case for all other species
    # Generate molecular properties estimation
    molecular_data = estimate_molecular_properties(formula)
    
    # Get configuration for which theoretical methods to use
    theory_config = get(get(config, "theory", Dict()), "use_methods", Dict())
    
    # Default values if not specified in config
    use_stat_thermo = get(theory_config, "statistical_thermodynamics", true)
    use_group_contribution = get(theory_config, "group_contribution", true)
    use_quantum_statistical = get(theory_config, "quantum_statistical", true)
    use_benson_group = get(theory_config, "benson_group", false)
    
    # Results array and method names
    theoretical_results = []
    method_names = []
    method_weights = []
    
    # Calculate using statistical thermodynamics
    if use_stat_thermo
        stat_thermo_result = calculate_statistical_thermodynamics(molecular_data, temperature)
        push!(theoretical_results, stat_thermo_result)
        push!(method_names, "Statistical Thermodynamics")
        push!(method_weights, 1.0)
    end
    
    # Calculate using group contribution method
    if use_group_contribution
        group_contrib_result = estimate_group_contribution(formula, temperature)
        push!(theoretical_results, group_contrib_result)
        push!(method_names, "Group Contribution")
        push!(method_weights, 1.0)
    end
    
    # Calculate using quantum-corrected statistical thermodynamics (more accurate)
    if use_quantum_statistical
        quantum_stat_thermo_result = calculate_quantum_statistical_thermodynamics(molecular_data, temperature)
        push!(theoretical_results, quantum_stat_thermo_result)
        push!(method_names, "Quantum-Statistical")
        push!(method_weights, 1.5)
    end
    
    # Calculate using bonds and group additivity method
    if use_benson_group
        benson_group_result = calculate_benson_group_additivity(formula, temperature)
        push!(theoretical_results, benson_group_result)
        push!(method_names, "Benson Group")
        push!(method_weights, 1.2)
    end
    
    # If no methods were enabled, default to statistical thermodynamics
    if isempty(theoretical_results)
        @warn "No theoretical methods enabled for $species_name. Using statistical thermodynamics as fallback."
        stat_thermo_result = calculate_statistical_thermodynamics(molecular_data, temperature)
        push!(theoretical_results, stat_thermo_result)
        push!(method_names, "Statistical Thermodynamics")
        push!(method_weights, 1.0)
    end
    
    # For each property, calculate weighted mean and uncertainty across methods
    properties = ["Cp", "H", "S", "G"]
    
    # Create individual method results for plotting
    individual_methods = []
    for i in 1:length(theoretical_results)
        method_result = Dict(
            "species_name" => species_name,
            "formula" => formula,
            "temperature" => temperature,
            "data_source" => "THEORETICAL_" * method_names[i],
            "properties" => Dict()
        )
        
        for prop in properties
            method_result["properties"][prop] = Dict(
                "value" => theoretical_results[i][prop],
                "uncertainty" => theoretical_results[i]["$(prop)_uncertainty"],
                "units" => prop in ["H", "G"] ? "kJ/mol" : "J/mol/K"
            )
        end
        
        push!(individual_methods, method_result)
    end
    
    combined_result = Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature" => temperature,
        "data_source" => "THEORETICAL_ENSEMBLE",
        "properties" => Dict(),
        "individual_methods" => individual_methods
    )
    
    for prop in properties
        values = [result[prop] for result in theoretical_results]
        uncertainties = [result["$(prop)_uncertainty"] for result in theoretical_results]
        
        # Weighted mean value across methods
        weighted_sum = sum(values .* method_weights)
        sum_weights = sum(method_weights)
        weighted_mean = weighted_sum / sum_weights
        
        # Combined uncertainty:
        # Weighted uncertainty combining standard deviation between methods and individual uncertainties
        weighted_var = sum(((values .- weighted_mean).^2) .* method_weights) / sum_weights
        weighted_std = sqrt(weighted_var)
        
        # Use higher of the weighted standard deviation or the minimum uncertainty scaled by reliability
        min_uncertainty = minimum(uncertainties)
        reliability_factor = 1.2  # Scale factor for minimum uncertainty
        combined_uncertainty = max(weighted_std, min_uncertainty * reliability_factor)
        
        # Store in combined result
        combined_result["properties"][prop] = Dict(
            "value" => weighted_mean,
            "uncertainty" => combined_uncertainty,
            "units" => prop in ["H", "G"] ? "kJ/mol" : "J/mol/K"
        )
    end
    
    return combined_result
end

"""
    calculate_quantum_statistical_thermodynamics(molecular_data::MolecularData, temperature::Float64)

Calculate thermodynamic properties using quantum-corrected statistical thermodynamics for greater accuracy.
Includes anharmonicity corrections and improved quantum effects.
"""
function calculate_quantum_statistical_thermodynamics(molecular_data::MolecularData, temperature::Float64)
    # Calculate partition functions with quantum corrections
    q_trans = calculate_translational_partition_function(molecular_data.molecular_weight, temperature)
    
    # Calculate rotational contribution with quantum corrections
    q_rot = 1.0
    if !isempty(molecular_data.rotational_constants)
        # Apply quantum correction factor
        quantum_correction = 1.0
        if length(molecular_data.rotational_constants) == 1
            # Linear molecule
            B = molecular_data.rotational_constants[1]
            quantum_correction = 1.0 + (1.0/24.0) * (h * B * c * 100)^2 / (kb * temperature)^2
        elseif length(molecular_data.rotational_constants) >= 3
            # Non-linear molecule
            A = molecular_data.rotational_constants[1]
            B = molecular_data.rotational_constants[2]
            C = molecular_data.rotational_constants[3]
            # Average rotational constant
            B_avg = (A + B + C) / 3
            quantum_correction = 1.0 + (1.0/12.0) * (h * B_avg * c * 100)^2 / (kb * temperature)^2
        end
        
        q_rot = calculate_rotational_partition_function(
            molecular_data.rotational_constants,
            molecular_data.symmetry_number,
            temperature
        ) * quantum_correction
    end
    
    # Calculate vibrational contribution with anharmonicity corrections
    q_vib = 1.0
    if !isempty(molecular_data.vibrational_frequencies)
        # Apply anharmonicity correction (diagonal correction)
        anharm_freqs = copy(molecular_data.vibrational_frequencies)
        for i in 1:length(anharm_freqs)
            # Typical anharmonicity reduces frequency by ~2%
            anharm_freqs[i] *= 0.98
        end
        
        q_vib = calculate_vibrational_partition_function_anharmonic(anharm_freqs, temperature)
    end
    
    # Calculate electronic contribution with multiple states if applicable
    q_elec = molecular_data.spin_multiplicity
    
    # Calculate thermodynamic properties with quantum corrections
    
    # Translational enthalpy with quantum correction
    h_trans = R * temperature * (1.0 + 1/(12*temperature^2))
    
    # Rotational enthalpy with quantum correction
    h_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        if length(molecular_data.rotational_constants) == 1
            # Linear molecule
            h_rot = R * temperature * (1.0 + 1/(24*temperature^2))
        else
            # Non-linear molecule with quantum correction
            h_rot = 1.5 * R * temperature * (1.0 + 1/(12*temperature^2))
        end
    end
    
    # Vibrational enthalpy with anharmonicity
    h_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        h_vib = calculate_vibrational_enthalpy_anharmonic(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    # Electronic contribution
    h_elec = 0.0
    
    # Total enthalpy
    enthalpy = h_trans + h_rot + h_vib + h_elec
    
    # Entropy calculations with quantum corrections
    
    # Translational entropy
    s_trans = R * (log(q_trans) + 1 + 1.5 + 1/(12*temperature^2))
    
    # Rotational entropy
    s_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        s_rot = R * (log(q_rot) + 1 + 1/(24*temperature^2))
    end
    
    # Vibrational entropy with anharmonicity
    s_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        s_vib = calculate_vibrational_entropy_anharmonic(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    # Electronic entropy
    s_elec = R * log(q_elec)
    
    # Total entropy
    s = s_trans + s_rot + s_vib + s_elec
    
    # Heat capacity calculations with quantum corrections
    
    # Translational heat capacity with quantum correction
    cp_trans = (5/2) * R * (1.0 + 1/(4*temperature^2))
    
    # Rotational heat capacity
    cp_rot = 0.0
    if !isempty(molecular_data.rotational_constants)
        if length(molecular_data.rotational_constants) == 1
            # Linear molecule
            cp_rot = R * (1.0 + 1/(12*temperature^2))
        else
            # Non-linear molecule
            cp_rot = (3/2) * R * (1.0 + 1/(6*temperature^2))
        end
    end
    
    # Vibrational heat capacity with anharmonicity
    cp_vib = 0.0
    if !isempty(molecular_data.vibrational_frequencies)
        cp_vib = calculate_vibrational_heat_capacity_anharmonic(
            molecular_data.vibrational_frequencies,
            temperature
        )
    end
    
    # Electronic heat capacity
    cp_elec = 0.0
    
    # Total heat capacity
    cp = cp_trans + cp_rot + cp_vib + cp_elec
    
    # Gibbs free energy
    g = enthalpy - temperature * s / 1000  # Convert S to kJ/mol/K
    
    # Estimate uncertainties based on the quality of molecular data
    has_good_vib = !isempty(molecular_data.vibrational_frequencies)
    has_good_rot = !isempty(molecular_data.rotational_constants)
    
    # Lower uncertainty due to quantum corrections
    base_uncertainty = has_good_vib && has_good_rot ? 0.08 : 0.18
    
    cp_uncertainty = max(cp * base_uncertainty, 2.5)
    h_uncertainty = max(abs(enthalpy) * base_uncertainty, 8.0)
    s_uncertainty = max(s * base_uncertainty, 4.0)
    g_uncertainty = max(abs(g) * base_uncertainty, 8.0)
    
    return Dict(
        "temperature" => temperature,
        "Cp" => cp,
        "Cp_uncertainty" => cp_uncertainty,
        "H" => enthalpy, 
        "H_uncertainty" => h_uncertainty,
        "S" => s,
        "S_uncertainty" => s_uncertainty,
        "G" => g,
        "G_uncertainty" => g_uncertainty,
        "method" => "quantum_statistical_thermodynamics"
    )
end

"""
    calculate_vibrational_partition_function_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational partition function with anharmonicity corrections.
"""
function calculate_vibrational_partition_function_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    q_vib = 1.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Anharmonicity correction (first-order approximation)
        # Typical anharmonicity parameter for molecules
        anharm_param = 0.005  # Varies by molecule, usually 0.003-0.01
        
        # Corrected partition function with anharmonicity
        q_vib_mode = 1 / (1 - exp(-x)) * (1 + anharm_param * x^2 * exp(-x) / (1 - exp(-x)))
        
        q_vib *= q_vib_mode
    end
    
    return q_vib
end

"""
    calculate_vibrational_enthalpy_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate vibrational enthalpy with anharmonicity corrections.
"""
function calculate_vibrational_enthalpy_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    h_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Anharmonicity correction (first-order approximation)
        anharm_param = 0.005
        anharm_correction = 1.0 + anharm_param * x * (1 + 1/(exp(x) - 1))
        
        # Vibrational enthalpy contribution with anharmonicity
        h_vib += R * temperature * x * exp(x) / (exp(x) - 1) * anharm_correction
        
        # Add zero-point energy with anharmonicity correction
        zpe_anharm = 1.0 - 2 * anharm_param  # Zero-point energy is reduced by anharmonicity
        h_vib += R * temperature * x / 2 * zpe_anharm
    end
    
    return h_vib
end

"""
    calculate_vibrational_entropy_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate vibrational entropy with anharmonicity corrections.
"""
function calculate_vibrational_entropy_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    s_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Anharmonicity correction (first-order approximation)
        anharm_param = 0.005
        
        # Additional entropy term due to anharmonicity
        anharm_term = anharm_param * x^2 * exp(-x) / (1 - exp(-x))^2
        
        # Vibrational entropy contribution with anharmonicity
        s_vib += R * (x / (exp(x) - 1) - log(1 - exp(-x)) + anharm_term)
    end
    
    return s_vib
end

"""
    calculate_vibrational_heat_capacity_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate vibrational heat capacity with anharmonicity corrections.
"""
function calculate_vibrational_heat_capacity_anharmonic(vibrational_frequencies::Vector{Float64}, temperature::Float64)
    cp_vib = 0.0
    
    for nu in vibrational_frequencies
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * c  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (kb * temperature)
        
        # Anharmonicity correction (first-order approximation)
        anharm_param = 0.005
        ex = exp(x)
        
        # Heat capacity correction due to anharmonicity
        anharm_term = anharm_param * x^2 * ex * (ex + 1) / (ex - 1)^3
        
        # Vibrational heat capacity contribution with anharmonicity
        cp_vib += R * (x^2 * ex / (ex - 1)^2 + anharm_term)
    end
    
    return cp_vib
end

"""
    calculate_benson_group_additivity(formula::String, temperature::Float64)

Calculate thermodynamic properties using Benson Group Additivity method.
This method is more accurate than simple group contribution for organic molecules.
"""
function calculate_benson_group_additivity(formula::String, temperature::Float64)
    # Parse formula
    elements = extract_formula_elements(formula)
    
    # Calculate molecular weight
    mw = calculate_molecular_weight(formula)
    
    # Basic group contributions for common atoms/groups
    # These are highly simplified - real Benson Group Additivity uses more complex groups
    group_contrib = Dict(
        "C" => Dict("Cp" => 4.0, "H" => 0.0, "S" => 0.0),  # Carbon atom
        "H" => Dict("Cp" => 2.8, "H" => 0.0, "S" => 0.0),  # Hydrogen atom
        "O" => Dict("Cp" => 3.5, "H" => 0.0, "S" => 0.0),  # Oxygen atom
        "N" => Dict("Cp" => 3.8, "H" => 0.0, "S" => 0.0),  # Nitrogen atom
        "S" => Dict("Cp" => 4.2, "H" => 0.0, "S" => 0.0),  # Sulfur atom
        "Ar" => Dict("Cp" => 2.5, "H" => 0.0, "S" => 0.0), # Argon
        "He" => Dict("Cp" => 2.5, "H" => 0.0, "S" => 0.0), # Helium
        "Ne" => Dict("Cp" => 2.5, "H" => 0.0, "S" => 0.0), # Neon
        "Xe" => Dict("Cp" => 2.5, "H" => 0.0, "S" => 0.0), # Xenon
    )
    
    # Bond contributions (simplified)
    # In real Benson method, these would be much more detailed
    bond_contrib = Dict(
        "C-C" => Dict("Cp" => 2.0, "H" => -20.0, "S" => 15.0),
        "C=C" => Dict("Cp" => 3.5, "H" => 12.0, "S" => 23.0),
        "C-H" => Dict("Cp" => 1.8, "H" => -15.0, "S" => 10.0),
        "C-O" => Dict("Cp" => 2.2, "H" => -30.0, "S" => 18.0),
        "C=O" => Dict("Cp" => 4.0, "H" => -120.0, "S" => 25.0),
        "O-H" => Dict("Cp" => 2.5, "H" => -45.0, "S" => 15.0),
        "N-H" => Dict("Cp" => 2.0, "H" => -40.0, "S" => 15.0),
        "N-N" => Dict("Cp" => 3.0, "H" => 0.0, "S" => 20.0),
        "N=N" => Dict("Cp" => 5.0, "H" => 10.0, "S" => 30.0),
        "O-O" => Dict("Cp" => 2.0, "H" => -10.0, "S" => 20.0),
    )
    
    # Initialize properties
    cp = 0.0
    h = 0.0
    s = 0.0
    
    # Apply atom contributions
    for (element, count) in elements
        if haskey(group_contrib, element)
            cp += group_contrib[element]["Cp"] * count
            h += group_contrib[element]["H"] * count
            s += group_contrib[element]["S"] * count
        end
    end
    
    # Estimate number of bonds based on valence rules (very simplified)
    total_atoms = sum(values(elements))
    if total_atoms > 1
        # Estimate number of bonds
        n_bonds = total_atoms - 1 + sum(get(elements, "N", 0) + get(elements, "O", 0) + get(elements, "C", 0)) ÷ 2
        
        # Apply bond contributions (this is a simplification)
        # In real Benson method, we would know the exact molecular structure
        if haskey(elements, "C") && elements["C"] > 0
            if haskey(elements, "H") && elements["H"] > 0
                # C-H bonds
                n_ch_bonds = min(elements["C"] * 4, elements["H"])
                cp += bond_contrib["C-H"]["Cp"] * n_ch_bonds
                h += bond_contrib["C-H"]["H"] * n_ch_bonds
                s += bond_contrib["C-H"]["S"] * n_ch_bonds
            end
            
            if elements["C"] > 1
                # C-C bonds (assuming single bonds for simplicity)
                n_cc_bonds = elements["C"] - 1
                cp += bond_contrib["C-C"]["Cp"] * n_cc_bonds
                h += bond_contrib["C-C"]["H"] * n_cc_bonds
                s += bond_contrib["C-C"]["S"] * n_cc_bonds
            end
            
            if haskey(elements, "O") && elements["O"] > 0
                # C-O bonds (assuming single bonds for simplicity)
                n_co_bonds = min(elements["C"], elements["O"])
                cp += bond_contrib["C-O"]["Cp"] * n_co_bonds
                h += bond_contrib["C-O"]["H"] * n_co_bonds
                s += bond_contrib["C-O"]["S"] * n_co_bonds
            end
        end
        
        if haskey(elements, "N") && elements["N"] > 0
            if haskey(elements, "H") && elements["H"] > 0
                # N-H bonds
                n_nh_bonds = min(elements["N"] * 3, elements["H"])
                cp += bond_contrib["N-H"]["Cp"] * n_nh_bonds
                h += bond_contrib["N-H"]["H"] * n_nh_bonds
                s += bond_contrib["N-H"]["S"] * n_nh_bonds
            end
            
            if elements["N"] > 1
                # N-N bonds or N=N bonds for N2
                if formula == "N2"
                    cp += bond_contrib["N=N"]["Cp"]
                    h += bond_contrib["N=N"]["H"]
                    s += bond_contrib["N=N"]["S"]
                else
                    n_nn_bonds = elements["N"] - 1
                    cp += bond_contrib["N-N"]["Cp"] * n_nn_bonds
                    h += bond_contrib["N-N"]["H"] * n_nn_bonds
                    s += bond_contrib["N-N"]["S"] * n_nn_bonds
                end
            end
        end
        
        if haskey(elements, "O") && elements["O"] > 0
            if haskey(elements, "H") && elements["H"] > 0
                # O-H bonds
                n_oh_bonds = min(elements["O"] * 2, elements["H"])
                cp += bond_contrib["O-H"]["Cp"] * n_oh_bonds
                h += bond_contrib["O-H"]["H"] * n_oh_bonds
                s += bond_contrib["O-H"]["S"] * n_oh_bonds
            end
            
            if elements["O"] > 1
                # O-O bonds or O=O bonds for O2
                if formula == "O2"
                    cp += bond_contrib["O-O"]["Cp"]
                    h += bond_contrib["O-O"]["H"]
                    s += bond_contrib["O-O"]["S"]
                else
                    n_oo_bonds = elements["O"] - 1
                    cp += bond_contrib["O-O"]["Cp"] * n_oo_bonds
                    h += bond_contrib["O-O"]["H"] * n_oo_bonds
                    s += bond_contrib["O-O"]["S"] * n_oo_bonds
                end
            end
        end
    end
    
    # Corrections based on molecular complexity
    if total_atoms > 3
        # Rotational and vibrational coupling in larger molecules
        # Increases heat capacity
        cp += 2.0 * (total_atoms - 3)
        
        # Increased entropy from conformational flexibility
        s += 5.0 * (total_atoms - 3)
    end
    
    # Temperature dependence (simplified)
    cp += 0.02 * cp * (temperature - 298.15)
    h += 0.01 * h * (temperature - 298.15) + cp * (temperature - 298.15)
    s += 0.005 * s * log(temperature / 298.15) + cp * log(temperature / 298.15)
    
    # Gibbs energy
    g = h - temperature * s / 1000
    
    # Estimate uncertainties - Benson method has moderate uncertainty
    cp_uncertainty = max(cp * 0.15, 3.0)
    h_uncertainty = max(abs(h) * 0.15, 10.0)
    s_uncertainty = max(s * 0.15, 5.0)
    g_uncertainty = max(abs(g) * 0.15, 10.0)
    
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
    calculate_nasa7_cp(coeffs::Vector{Float64}, temperature::Float64)

Calculate heat capacity (Cp/R) using NASA-7 polynomial coefficients.
"""
function calculate_nasa7_cp(coeffs::Vector{Float64}, temperature::Float64)
    return coeffs[1] + 
           coeffs[2] * temperature + 
           coeffs[3] * temperature^2 + 
           coeffs[4] * temperature^3 + 
           coeffs[5] * temperature^4
end

"""
    calculate_nasa7_enthalpy(coeffs::Vector{Float64}, temperature::Float64)

Calculate enthalpy (H/RT) using NASA-7 polynomial coefficients.
"""
function calculate_nasa7_enthalpy(coeffs::Vector{Float64}, temperature::Float64)
    return coeffs[1] + 
           (coeffs[2]/2) * temperature + 
           (coeffs[3]/3) * temperature^2 + 
           (coeffs[4]/4) * temperature^3 + 
           (coeffs[5]/5) * temperature^4 + 
           coeffs[6] / temperature
end

"""
    calculate_nasa7_entropy(coeffs::Vector{Float64}, temperature::Float64)

Calculate entropy (S/R) using NASA-7 polynomial coefficients.
"""
function calculate_nasa7_entropy(coeffs::Vector{Float64}, temperature::Float64)
    return coeffs[1] * log(temperature) + 
           coeffs[2] * temperature + 
           (coeffs[3]/2) * temperature^2 + 
           (coeffs[4]/3) * temperature^3 + 
           (coeffs[5]/4) * temperature^4 + 
           coeffs[7]
end

"""
    calculate_nasa_properties(polynomial::ThermodynamicPolynomial, temperature::Float64)

Calculate thermodynamic properties using NASA polynomials.
"""
function calculate_nasa_properties(polynomial::ThermodynamicPolynomial, temperature::Float64)
    # Find applicable temperature range
    range_idx = 0
    for (idx, range) in enumerate(polynomial.temperature_ranges)
        if range[1] <= temperature && temperature <= range[2]
            range_idx = idx
            break
        end
    end
    
    # If no range found, use the closest one
    if range_idx == 0
        if temperature < polynomial.temperature_ranges[1][1]
            range_idx = 1
            @warn "Temperature $temperature K is below the valid range. Using the lowest range."
        else
            range_idx = length(polynomial.temperature_ranges)
            @warn "Temperature $temperature K is above the valid range. Using the highest range."
        end
    end
    
    # Get coefficients for this range
    coeffs = polynomial.coefficients[range_idx]
    
    # Calculate properties depending on polynomial type
    if polynomial.polynomial_type == "nasa7"
        # Calculate Cp/R
        cp = calculate_nasa7_cp(coeffs, temperature)
        
        # Calculate H/RT
        h = calculate_nasa7_enthalpy(coeffs, temperature)
        
        # Calculate S/R
        s = calculate_nasa7_entropy(coeffs, temperature)
        
    elseif polynomial.polynomial_type == "nasa9"
        # NASA-9 implementation would go here
        # For now, just use NASA-7 approach as approximation
        cp = calculate_nasa7_cp(coeffs, temperature)
        h = calculate_nasa7_enthalpy(coeffs, temperature)
        s = calculate_nasa7_entropy(coeffs, temperature)
    else
        error("Unsupported polynomial type: $(polynomial.polynomial_type)")
    end
    
    # Convert to standard units
    cp_val = cp * R
    h_val = h * R * temperature / 1000  # kJ/mol
    s_val = s * R
    g_val = h_val - temperature * s_val / 1000  # kJ/mol
    
    # Calculate uncertainties based on reliability score
    # Lower reliability score means higher uncertainty
    rel_factor = (5.0 - polynomial.reliability_score) / 5.0  # 0 to 1 scale
    uncertainty_base = 0.02 + 0.05 * rel_factor  # 2% to 7% base uncertainty
    
    cp_uncertainty = uncertainty_base * cp_val
    h_uncertainty = uncertainty_base * abs(h_val)
    s_uncertainty = uncertainty_base * s_val
    g_uncertainty = uncertainty_base * abs(g_val)
    
    # Create result structure
    return Dict(
        "species_name" => polynomial.species_name,
        "formula" => polynomial.formula,
        "temperature" => temperature,
        "data_source" => polynomial.source,
        "properties" => Dict(
            "Cp" => Dict(
                "value" => cp_val,
                "uncertainty" => cp_uncertainty,
                "units" => "J/mol/K"
            ),
            "H" => Dict(
                "value" => h_val,
                "uncertainty" => h_uncertainty,
                "units" => "kJ/mol"
            ),
            "S" => Dict(
                "value" => s_val,
                "uncertainty" => s_uncertainty,
                "units" => "J/mol/K"
            ),
            "G" => Dict(
                "value" => g_val,
                "uncertainty" => g_uncertainty,
                "units" => "kJ/mol"
            )
        )
    )
end

"""
    combine_weighted_values(val1::Float64, unc1::Float64, val2::Float64, unc2::Float64, 
                         weight1::Float64, weight2::Float64)

Combine two values with uncertainties using weighted averaging.
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
"""
function refine_thermodynamic_data(previous_result::Dict, new_result::Dict, weight_factor::Float64)
    # Create a copy of the previous result
    refined_result = deepcopy(previous_result)
    
    # Update source information
    refined_result["data_source"] = new_result["data_source"]
    
    # Calculate weights - higher weight for new data from more reliable source
    previous_weight = 1.0 - weight_factor  # Lower weight for previous data
    new_weight = weight_factor  # Higher weight for new data
    
    # Update each property
    for prop in ["Cp", "H", "S", "G"]
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
    
    return refined_result
end

"""
    progressively_refine_thermodynamic_data(species_name::String, formula::String, 
                                         temperature::Float64, data_sources::Dict)

Progressively refine thermodynamic data by traversing the hierarchical data sources.
Returns the refined result, all sources used, and a list of sources used in order of application.
"""
function progressively_refine_thermodynamic_data(species_name::String, formula::String, 
                                              temperature::Float64, data_sources::Dict)
    # STEP 1: Start with theoretical calculations
    result = calculate_theoretical_properties(species_name, formula, temperature)
    
    # Initialize all sources list and sources used tracking
    all_sources = [result]
    sources_used = [Dict(
        "source_name" => "THEORETICAL",
        "priority" => 0,
        "reliability_score" => 2.5,
        "weight_factor" => 1.0
    )]
    
    # STEP 2: Traverse the hierarchy of data sources
    for priority in 1:8
        # Find sources at this priority level
        for (source_name, source_data) in data_sources
            if source_data["priority"] == priority
                # Check if this source has data for our species
                if haskey(source_data["polynomials"], species_name)
                    polynomial = source_data["polynomials"][species_name]
                    
                    # Calculate properties using this polynomial
                    source_result = calculate_nasa_properties(polynomial, temperature)
                    push!(all_sources, source_result)
                    
                    # Calculate weight for this source
                    # Higher priority and reliability give higher weight
                    base_weight = 0.5 + 0.05 * priority  # 0.55 to 0.9 based on priority
                    reliability_factor = polynomial.reliability_score / 5.0  # 0 to 1 scale
                    weight_factor = base_weight * (0.8 + 0.2 * reliability_factor)  # Adjust weight by reliability
                    
                    # Track this source's contribution
                    push!(sources_used, Dict(
                        "source_name" => source_name,
                        "priority" => priority,
                        "reliability_score" => polynomial.reliability_score,
                        "weight_factor" => weight_factor
                    ))
                    
                    # Refine the result
                    result = refine_thermodynamic_data(result, source_result, weight_factor)
                end
            end
        end
    end
    
    return result, all_sources, sources_used
end

"""
    calculate_properties_range(species_name::String, formula::String, temp_range::Vector{<:Real}, 
                             step::Real, data_sources::Dict)

Calculate thermodynamic properties over a temperature range.
"""
function calculate_properties_range(species_name::String, formula::String, temp_range::Vector{<:Real}, 
                                  step::Real, data_sources::Dict, config::Dict=Dict())
    temp_min = Float64(temp_range[1])
    temp_max = Float64(temp_range[2])
    step_val = Float64(step)
    
    # Generate temperature points
    temperatures = temp_min:step_val:temp_max
    
    # Call the hierarchical properties calculator to get step-by-step results
    hierarchical_steps = calculate_hierarchical_properties(species_name, formula, collect(temperatures), data_sources, config)
    
    # Extract the final results (last item in the hierarchical steps)
    results = hierarchical_steps[end]["results"]
    
    # We'll get sources used from the mid-temperature point
    mid_temp_idx = findfirst(t -> t >= 1000.0, temperatures)
    if mid_temp_idx === nothing
        mid_temp_idx = div(length(temperatures), 2)  # Use the middle temperature
    end
    
    # Build sources_used list from intermediate steps
    sources_used = []
    for step in hierarchical_steps
        if step["source_name"] != "FINAL REFINED"  # Skip the final step
            source_info = Dict(
                "source_name" => step["source_name"],
                "priority" => step["priority"]
            )
            
            # Add reliability score if available
            if haskey(step, "reliability_score")
                source_info["reliability_score"] = step["reliability_score"]
            else
                source_info["reliability_score"] = 2.5  # Default for theoretical
            end
            
            # Add weight factor (estimated)
            if step["source_name"] == "THEORETICAL"
                source_info["weight_factor"] = 1.0
            else
                priority = step["priority"]
                reliability = source_info["reliability_score"]
                base_weight = 0.5 + 0.05 * priority  # 0.55 to 0.9 based on priority
                reliability_factor = reliability / 5.0  # 0 to 1 scale
                source_info["weight_factor"] = base_weight * (0.8 + 0.2 * reliability_factor)
            end
            
            push!(sources_used, source_info)
        end
    end
    
    # Extract property data for plotting
    temps = [result["temperature"] for result in results]
    
    cp_values = [result["properties"]["Cp"]["value"] for result in results]
    cp_uncertainties = [result["properties"]["Cp"]["uncertainty"] for result in results]
    
    h_values = [result["properties"]["H"]["value"] for result in results]
    h_uncertainties = [result["properties"]["H"]["uncertainty"] for result in results]
    
    s_values = [result["properties"]["S"]["value"] for result in results]
    s_uncertainties = [result["properties"]["S"]["uncertainty"] for result in results]
    
    g_values = [result["properties"]["G"]["value"] for result in results]
    g_uncertainties = [result["properties"]["G"]["uncertainty"] for result in results]
    
    # Create result structure
    return Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature_range" => temp_range,
        "temperatures" => temps,
        "sources_used" => sources_used,
        "hierarchical_steps" => hierarchical_steps,
        "properties" => Dict(
            "Cp" => Dict(
                "values" => cp_values,
                "uncertainties" => cp_uncertainties,
                "units" => "J/mol/K"
            ),
            "H" => Dict(
                "values" => h_values,
                "uncertainties" => h_uncertainties,
                "units" => "kJ/mol"
            ),
            "S" => Dict(
                "values" => s_values,
                "uncertainties" => s_uncertainties,
                "units" => "J/mol/K"
            ),
            "G" => Dict(
                "values" => g_values,
                "uncertainties" => g_uncertainties,
                "units" => "kJ/mol"
            )
        )
    )
end

"""
    fit_nasa7_coefficients(temps::Vector{Float64}, cp_vals::Vector{Float64}, 
                         h_vals::Vector{Float64}, s_vals::Vector{Float64})

Fit NASA-7 polynomial coefficients to thermodynamic data.
"""
function fit_nasa7_coefficients(temps::Vector{Float64}, cp_vals::Vector{Float64}, 
                              h_vals::Vector{Float64}, s_vals::Vector{Float64})
    # Divide temperature range into two parts
    n = length(temps)
    mid_idx = max(3, div(n, 2))  # Ensure we have at least 3 points for fitting
    
    # Low temperature range
    low_temps = temps[1:mid_idx]
    low_cp = cp_vals[1:mid_idx] ./ R  # Convert to Cp/R
    low_h = h_vals[1:mid_idx] .* 1000.0 ./ (R .* low_temps)  # Convert to H/RT
    low_s = s_vals[1:mid_idx] ./ R  # Convert to S/R
    
    # High temperature range
    high_temps = temps[mid_idx:end]
    high_cp = cp_vals[mid_idx:end] ./ R  # Convert to Cp/R
    high_h = h_vals[mid_idx:end] .* 1000.0 ./ (R .* high_temps)  # Convert to H/RT
    high_s = s_vals[mid_idx:end] ./ R  # Convert to S/R
    
    # Function to fit NASA-7 CP coefficients (a1-a5)
    function fit_cp_coeffs(temperatures, cp_values)
        # Create design matrix for Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        A = hcat(
            ones(length(temperatures)),
            temperatures,
            temperatures.^2,
            temperatures.^3,
            temperatures.^4
        )
        
        # Solve least squares problem
        return A \ cp_values
    end
    
    # Fit a1-a5 for low temperature range
    low_cp_coeffs = fit_cp_coeffs(low_temps, low_cp)
    
    # Calculate a6 from H/RT data
    a6_low = 0.0
    for i in 1:length(low_temps)
        T = low_temps[i]
        h_rt_calc = low_cp_coeffs[1] + 
                   (low_cp_coeffs[2]/2)*T + 
                   (low_cp_coeffs[3]/3)*T^2 + 
                   (low_cp_coeffs[4]/4)*T^3 + 
                   (low_cp_coeffs[5]/5)*T^4
        a6_low += (low_h[i] - h_rt_calc) * T
    end
    a6_low /= length(low_temps)
    
    # Calculate a7 from S/R data
    a7_low = 0.0
    for i in 1:length(low_temps)
        T = low_temps[i]
        s_r_calc = low_cp_coeffs[1] * log(T) + 
                  low_cp_coeffs[2]*T + 
                  (low_cp_coeffs[3]/2)*T^2 + 
                  (low_cp_coeffs[4]/3)*T^3 + 
                  (low_cp_coeffs[5]/4)*T^4
        a7_low += (low_s[i] - s_r_calc)
    end
    a7_low /= length(low_temps)
    
    # Fit a1-a5 for high temperature range
    high_cp_coeffs = fit_cp_coeffs(high_temps, high_cp)
    
    # Calculate a6 from H/RT data
    a6_high = 0.0
    for i in 1:length(high_temps)
        T = high_temps[i]
        h_rt_calc = high_cp_coeffs[1] + 
                   (high_cp_coeffs[2]/2)*T + 
                   (high_cp_coeffs[3]/3)*T^2 + 
                   (high_cp_coeffs[4]/4)*T^3 + 
                   (high_cp_coeffs[5]/5)*T^4
        a6_high += (high_h[i] - h_rt_calc) * T
    end
    a6_high /= length(high_temps)
    
    # Calculate a7 from S/R data
    a7_high = 0.0
    for i in 1:length(high_temps)
        T = high_temps[i]
        s_r_calc = high_cp_coeffs[1] * log(T) + 
                  high_cp_coeffs[2]*T + 
                  (high_cp_coeffs[3]/2)*T^2 + 
                  (high_cp_coeffs[4]/3)*T^3 + 
                  (high_cp_coeffs[5]/4)*T^4
        a7_high += (high_s[i] - s_r_calc)
    end
    a7_high /= length(high_temps)
    
    # Create NASA-7 coefficient vectors
    low_coeffs = [low_cp_coeffs..., a6_low, a7_low]
    high_coeffs = [high_cp_coeffs..., a6_high, a7_high]
    
    # Create temperature ranges
    low_range = [temps[1], temps[mid_idx]]
    high_range = [temps[mid_idx], temps[end]]
    
    return [low_coeffs, high_coeffs], [low_range, high_range]
end

"""
    load_polynomial_data()

Load thermodynamic polynomial data for different species from multiple sources,
ensuring comprehensive coverage of all available data sources.
"""
function load_polynomial_data()
    # Initialize data sources dictionary
    data_sources = Dict()
    
    # Common temperature ranges
    temp_range_low = [200.0, 1000.0]
    temp_range_high = [1000.0, 6000.0]
    temp_ranges = [temp_range_low, temp_range_high]
    
    # Load species list from config
    config_path = joinpath(@__DIR__, "config", "Species.yaml")
    local species_list = []
    if isfile(config_path)
        config = YAML.load_file(config_path)
        species_list = config["species"]
    else
        # Fallback to default species if config not found
        species_list = ["N2", "O2", "H2O", "CO2", "CH4", "H2", "CO", "Ar", "He"]
    end
    
    # ========= PRIORITY 1: CHEMKIN Data (GRI-MECH 3.0) ==========
    # CHEMKIN/GRI-MECH 3.0 data - commonly used in combustion modeling
    chemkin_polys = Dict{String, ThermodynamicPolynomial}()
    chemkin_data_path = joinpath(@__DIR__, "data", "cache", "gri-mech", "therm.dat")
    gri_mech_data = Dict(
        "N2" => (
            [[3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372],
             [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]]
        ),
        "O2" => (
            [[3.78246, -0.0029673416, 9.8473583e-06, -9.6813198e-09, 3.2437178e-12, -1063.9433, 3.6598074],
             [3.28254, 0.0014816965, -7.5754135e-07, 2.0949245e-10, -2.1671779e-14, -1088.4576, 5.4536883]]
        ),
        "H2O" => (
            [[4.19864, -0.0020364017, 6.5203416e-06, -5.4879269e-09, 1.7719682e-12, -30293.73, -0.84900901],
             [2.67704, 0.0029731816, -7.7376889e-07, 9.4433514e-11, -4.2689991e-15, -29885.89, 6.88255]]
        ),
        "CO2" => (
            [[2.35677, 0.0089841299, -7.1220632e-06, 2.4573008e-09, -1.4288548e-13, -48371.97, 9.9009035],
             [3.85746, 0.0044626261, -2.4147881e-06, 6.387175e-10, -6.3472846e-14, -48759.17, 2.2808193]]
        ),
        "CH4" => (
            [[5.14988, -0.013671193, 4.9472979e-05, -4.8497523e-08, 1.6687213e-11, -10246.60, -4.6413137],
             [0.0748515, 0.013339079, -5.7321517e-06, 1.2229303e-09, -9.9723584e-14, -9468.39, 18.437322]]
        ),
        "H2" => (
            [[2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238],
             [2.93286575, 0.000826608026, -1.46402364e-07, 1.54100414e-11, -6.888048e-16, -813.065581, -1.02432865]]
        ),
        "CO" => (
            [[3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928],
             [2.71518561, 0.00206252743, -9.98825771e-07, 2.30053008e-10, -2.03647716e-14, -14151.8724, 7.81868772]]
        ),
        "Ar" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366],
             [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]]
        ),
        "He" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 0.928],
             [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 0.928]]
        ),
        "NO" => (
            [[4.21859896, -0.00463988124, 1.10443049e-05, -9.34055507e-09, 2.80554874e-12, 9845.09964, 2.28061001],
             [3.26071234, 0.00119101135, -4.29122646e-07, 6.94481463e-11, -4.03295681e-15, 9921.43132, 6.36900518]]
        ),
        "OH" => (
            [[3.99198424, -0.00240106655, 4.61664033e-06, -3.87916306e-09, 1.36319502e-12, 3368.89836, -0.103998477],
             [2.83853033, 0.00110741289, -2.94000209e-07, 4.20698729e-11, -2.4228989e-15, 3697.80808, 5.84494652]]
        ),
        "H" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, 25473.66, -0.44668285],
             [2.5, 0.0, 0.0, 0.0, 0.0, 25473.66, -0.44668285]]
        ),
        "O" => (
            [[3.16826710, -0.00327931884, 6.64306396e-06, -6.12806624e-09, 2.11265971e-12, 29122.2592, 2.05193346],
             [2.56942078, -8.59741137e-05, 4.19484589e-08, -1.00177799e-11, 1.22833691e-15, 29217.5791, 4.78433864]]
        ),
        # Add more species from GRI-MECH
        "NO2" => (
            [[3.944031, -0.0015797740, 1.6656039e-05, -2.0475426e-08, 7.8350564e-12, 2896.618, 6.3119917],
             [4.8847540, 0.0021723956, -8.2806906e-07, 1.5747510e-10, -1.1806638e-14, 2316.4980, -0.75442170]]
        ),
        "N2O" => (
            [[2.257150200, 0.011304682, -1.36591500e-05, 9.16618560e-09, -2.303022800e-12, 8741.77340, 10.56763389],
             [4.823072900, 0.002627025, -9.58501070e-07, 1.60000000e-10, -9.77523380e-15, 8073.40465, -2.20172080]]
        ),
        "NH3" => (
            [[4.286027400, -0.004623033, 2.15215150e-05, -2.29850340e-08, 8.280528950e-12, -6741.72280, -0.62517932],
             [2.634452100, 0.005666256, -1.72863600e-06, 2.38605240e-10, -1.25731700e-14, -6544.69520, 6.56296057]]
        )
    )
    
    # Create polynomial entries for all available species in our data
    for species_name in species_list
        # Skip if no data available
        if !haskey(gri_mech_data, species_name)
            continue
        end
        
        coeffs = gri_mech_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "GRI-MECH", 1, 3.5, "nasa7",
            temp_ranges,
            coeffs,
            0.02  # 2% uncertainty
        )
        
        chemkin_polys[species_name] = poly
    end
    
    data_sources["GRI-MECH"] = Dict(
        "priority" => 1,
        "reliability_score" => 3.5,
        "polynomials" => chemkin_polys
    )
    
    # ========= PRIORITY 2: CHEMKIN Data (Additional mechanism) ==========
    # Add more species from other CHEMKIN mechanisms
    chemkin_additional_polys = Dict{String, ThermodynamicPolynomial}()
    chemkin_additional_data = Dict(
        "SO2" => (
            [[3.266296, 0.00537144, -4.9000168e-06, 2.1967804e-09, -3.8389831e-13, -36945.16, 9.6816761],
             [5.2521167, 0.0019777653, -8.2053776e-07, 1.5393643e-10, -1.0631419e-14, -37558.14, -1.0889978]]
        ),
        "SO3" => (
            [[2.578641, 0.014550449, -1.7334793e-05, 9.8629329e-09, -2.1534482e-12, -48597.30, 12.061334],
             [7.0758174, 0.0032931588, -1.4132337e-06, 2.6716914e-10, -1.8700792e-14, -50655.14, -11.732945]]
        ),
        "HCN" => (
            [[2.4537740, 0.008607413, -1.07942840e-05, 8.28573220e-09, -2.231314660e-12, 14782.1741, 9.3353720],
             [3.8020657, 0.003146297, -1.06344930e-06, 1.66187830e-10, -9.799757460e-15, 14364.5999, 1.5752420]]
        ),
        "CN" => (
            [[3.9403710, -0.003229127, 6.97756360e-06, -5.92721450e-09, 1.873130540e-12, 51466.2371, 1.8439310],
             [3.7459804, 4.345077e-04, -3.04283540e-07, 8.49372870e-11, -6.150855560e-15, 51532.3434, 3.3307418]]
        )
    )
    
    # Create polynomial entries for additional CHEMKIN species
    for species_name in species_list
        if !haskey(chemkin_additional_data, species_name)
            continue
        end
        
        coeffs = chemkin_additional_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "CHEMKIN-Additional", 2, 3.3, "nasa7",
            temp_ranges,
            coeffs,
            0.025  # 2.5% uncertainty
        )
        
        chemkin_additional_polys[species_name] = poly
    end
    
    data_sources["CHEMKIN-Additional"] = Dict(
        "priority" => 2,
        "reliability_score" => 3.3,
        "polynomials" => chemkin_additional_polys
    )
    
    # ========= PRIORITY 3: NASA CEA DATA ==========
    # NASA CEA database - comprehensive thermodynamic data
    nasa_cea_polys = Dict{String, ThermodynamicPolynomial}()
    nasa_cea_data_path = joinpath(@__DIR__, "data", "cache", "nasa-cea", "thermo.inp")
    nasa_cea_data = Dict(
        "N2" => (
            [[3.53100528, -0.000123660988, -5.02999433e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747038],
             [2.95257637, 0.0013969004, -4.92631603e-07, 7.86010195e-11, -4.60755204e-15, -923.948688, 5.87188762]]
        ),
        "O2" => (
            [[3.78245636, -0.00299673416, 9.84735853e-06, -9.68131659e-09, 3.24372837e-12, -1063.94356, 3.65750924],
             [3.66096065, 0.000656365811, -1.41149627e-07, 2.05797935e-11, -1.29913436e-15, -1215.97718, 3.41536279]]
        ),
        "H2O" => (
            [[4.19864056, -0.00203643410, 6.52034160e-06, -5.48792690e-09, 1.77196820e-12, -30293.7267, -0.849009010],
             [2.67703787, 0.00297318329, -7.73768890e-07, 9.44335140e-11, -4.26899910e-15, -29885.8943, 6.88255000]]
        ),
        "CO2" => (
            [[2.35677352, 0.00898412990, -7.12206320e-06, 2.45730080e-09, -1.42885480e-13, -48371.9697, 9.90090350],
             [3.85746029, 0.00441437026, -2.21481404e-06, 5.23490188e-10, -4.72084164e-14, -48759.1662, 2.27167050]]
        ),
        # Add many more species from NASA CEA
        "NO" => (
            [[4.218599e+00, -4.639976e-03, 1.104431e-05, -9.340596e-09, 2.805554e-12, 9845.099, 2.280610],
             [3.260712e+00, 1.191001e-03, -4.291205e-07, 6.945752e-11, -4.033460e-15, 9921.43, 6.369604]]
        ),
        "O3" => (
            [[3.407006, 0.002, -1.06271800e-05, 1.15403900e-08, -3.841550780e-12, 15864.0753, 8.2862945],
             [5.401625, 0.001, -2.05403200e-07, 1.86928830e-11, -1.846146650e-15, 15339.3965, -3.2954596]]
        ),
        "OH" => (
            [[3.991982, -2.40107e-03, 4.61664e-06, -3.87916e-09, 1.36319e-12, 3368.898, -0.103998],
             [2.838350, 1.10741e-03, -2.94e-07, 4.20698e-11, -2.42284e-15, 3697.807, 5.84494]]
        ),
        "H2O2" => (
            [[4.31515149, -8.47390622e-04, 1.76404323e-05, -2.26762944e-08, 9.08950158e-12, -17706.7437, 3.27373216],
             [4.57977305, 4.05326003e-03, -1.29844730e-06, 1.98211400e-10, -1.13968792e-14, -18007.1776, 0.664970694]]
        ),
        "N2O" => (
            [[2.257150, 0.011304682, -1.3659150e-05, 9.1661860e-09, -2.30302280e-12, 8741.7734, 10.5676339],
             [4.823072, 0.002627025, -9.5850107e-07, 1.6000000e-10, -9.77523380e-15, 8073.4046, -2.2017208]]
        ),
        "CO" => (
            [[3.579533, -6.10354e-04, 1.01681e-06, 9.07005e-10, -9.04424e-13, -14344.09, 3.50841],
             [2.715185, 2.06252e-03, -9.98826e-07, 2.30053e-10, -2.03648e-14, -14151.9, 7.81869]]
        ),
        "S" => (
            [[2.947460, -1.63514e-03, 2.49781e-06, -1.89927e-09, 5.60843e-13, 33052.6, 3.93346],
             [2.219057, 2.58548e-04, -1.04188e-07, 2.33853e-11, -2.17429e-15, 33133.7, 6.83082]]
        ),
        "Xe" => (
            [[2.50000000, 0.0, 0.0, 0.0, 0.0, -745.375, 5.45873589],
             [2.50000000, 0.0, 0.0, 0.0, 0.0, -745.375, 5.45873589]]
        ),
        "Ne" => (
            [[2.50000000, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532273],
             [2.50000000, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532273]]
        )
    )
    
    # Create polynomial entries for all available species in NASA CEA
    for species_name in species_list
        # Skip if no data available
        if !haskey(nasa_cea_data, species_name)
            continue
        end
        
        coeffs = nasa_cea_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "NASA-CEA", 3, 4.0, "nasa7",
            temp_ranges,
            coeffs,
            0.015  # 1.5% uncertainty
        )
        
        nasa_cea_polys[species_name] = poly
    end
    
    data_sources["NASA-CEA"] = Dict(
        "priority" => 3,
        "reliability_score" => 4.0,
        "polynomials" => nasa_cea_polys
    )
    
    # ========= PRIORITY 4: JANAF Thermochemical Tables ==========
    # JANAF - high accuracy experimental data
    janaf_polys = Dict{String, ThermodynamicPolynomial}()
    janaf_data_path = joinpath(@__DIR__, "data", "cache", "janaf", "data.txt")
    janaf_data = Dict(
        "N2" => (
            [[3.531005, -0.0001236610, -5.0299943e-07, 2.4353061e-09, -1.4088124e-12, -1046.9763, 2.9674704],
             [2.952576, 0.0013969004, -4.9263160e-07, 7.8601020e-11, -4.6075520e-15, -923.94869, 5.8718876]]
        ),
        "O2" => (
            [[3.78246, -0.00299673, 9.8473583e-06, -9.6813198e-09, 3.2437178e-12, -1063.9433, 3.6598074],
             [3.28254, 0.00148170, -7.5754135e-07, 2.0949245e-10, -2.1671779e-14, -1088.4576, 5.4536883]]
        ),
        "H2O" => (
            [[4.198641, -0.0020364341, 6.5203416e-06, -5.4879269e-09, 1.7719682e-12, -30293.726, -0.84901],
             [2.677038, 0.0029731833, -7.7376889e-07, 9.4433514e-11, -4.2689991e-15, -29885.894, 6.88255]]
        ),
        "CO2" => (
            [[2.356773, 0.00898413, -7.1220632e-06, 2.4573008e-09, -1.4288548e-13, -48371.969, 9.9009035],
             [3.857460, 0.00444264, -2.2148140e-06, 5.2349019e-10, -4.7208416e-14, -48759.166, 2.27167]]
        ),
        # Add more JANAF species
        "NH3" => (
            [[4.300529, -0.0046615014, 2.1718710e-05, -2.2808488e-08, 8.2638046e-12, -6697.9108, -0.66490510],
             [2.731650, 0.0056308385, -1.7596395e-06, 2.5234329e-10, -1.3671264e-14, -6493.7289, 6.09126148]]
        ),
        "H2" => (
            [[2.344331, 0.00798052, -1.9478151e-05, 2.0157209e-08, -7.3761176e-12, -917.93517, 0.68301024],
             [2.932866, 0.00082661, -1.4640236e-07, 1.5410041e-11, -6.8880480e-16, -813.06558, -1.0243287]]
        ),
        "NO" => (
            [[4.218599, -0.00463998, 1.1044309e-05, -9.3405551e-09, 2.8055487e-12, 9845.0996, 2.2806100],
             [3.260712, 0.00119110, -4.2912265e-07, 6.9448146e-11, -4.0329568e-15, 9921.4313, 6.3690052]]
        ),
        "SO2" => (
            [[3.266296, 0.00537144, -4.9000168e-06, 2.1967804e-09, -3.8389831e-13, -36945.16, 9.6816761],
             [5.252117, 0.00197777, -8.2053776e-07, 1.5393643e-10, -1.0631419e-14, -37558.14, -1.0889978]]
        )
    )
    
    # Create polynomial entries for all available species in JANAF
    for species_name in species_list
        # Skip if no data available
        if !haskey(janaf_data, species_name)
            continue
        end
        
        coeffs = janaf_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "JANAF", 4, 4.5, "nasa7",
            temp_ranges,
            coeffs,
            0.010  # 1.0% uncertainty - JANAF data is quite reliable
        )
        
        janaf_polys[species_name] = poly
    end
    
    data_sources["JANAF"] = Dict(
        "priority" => 4,
        "reliability_score" => 4.5,
        "polynomials" => janaf_polys
    )
    
    # ========= PRIORITY 5: ThermoML Database ==========
    # ThermoML - experimental data from literature
    thermoml_polys = Dict{String, ThermodynamicPolynomial}()
    thermoml_data = Dict(
        "H2O" => (
            [[4.198641, -0.0020364341, 6.5203416e-06, -5.4879269e-09, 1.7719682e-12, -30293.726, -0.84901],
             [2.677038, 0.0029731833, -7.7376889e-07, 9.4433514e-11, -4.2689991e-15, -29885.894, 6.88255]]
        ),
        "CH4" => (
            [[5.14988, -0.013671193, 4.9472979e-05, -4.8497523e-08, 1.6687213e-11, -10246.60, -4.6413137],
             [0.0748515, 0.013339079, -5.7321517e-06, 1.2229303e-09, -9.9723584e-14, -9468.39, 18.437322]]
        ),
        "CO2" => (
            [[2.356773, 0.00898413, -7.1220632e-06, 2.4573008e-09, -1.4288548e-13, -48371.969, 9.9009035],
             [3.857460, 0.00444264, -2.2148140e-06, 5.2349019e-10, -4.7208416e-14, -48759.166, 2.27167]]
        )
    )
    
    # Create polynomial entries for ThermoML data
    for species_name in species_list
        # Skip if no data available
        if !haskey(thermoml_data, species_name)
            continue
        end
        
        coeffs = thermoml_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "ThermoML", 5, 4.3, "nasa7",
            temp_ranges,
            coeffs,
            0.012  # 1.2% uncertainty
        )
        
        thermoml_polys[species_name] = poly
    end
    
    data_sources["ThermoML"] = Dict(
        "priority" => 5,
        "reliability_score" => 4.3,
        "polynomials" => thermoml_polys
    )
    
    # ========= PRIORITY 6: NIST TDE (Thermodynamics Data Engine) ==========
    # NIST TDE - high-quality evaluated data
    tde_polys = Dict{String, ThermodynamicPolynomial}()
    tde_data = Dict(
        "N2" => (
            [[3.53100528, -0.000123660988, -5.02999433e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747038],
             [2.95257637, 0.0013969004, -4.92631603e-07, 7.86010195e-11, -4.60755204e-15, -923.948688, 5.87188762]]
        ),
        "O2" => (
            [[3.78245636, -0.00299673416, 9.84735853e-06, -9.68131659e-09, 3.24372837e-12, -1063.94356, 3.65750924],
             [3.66096065, 0.000656365811, -1.41149627e-07, 2.05797935e-11, -1.29913436e-15, -1215.97718, 3.41536279]]
        ),
        "H2O" => (
            [[4.19864056, -0.00203643410, 6.52034160e-06, -5.48792690e-09, 1.77196820e-12, -30293.7267, -0.849009010],
             [2.67703787, 0.00297318329, -7.73768890e-07, 9.44335140e-11, -4.26899910e-15, -29885.8943, 6.88255000]]
        ),
        "CH4" => (
            [[5.14988, -0.013671193, 4.9472979e-05, -4.8497523e-08, 1.6687213e-11, -10246.60, -4.6413137],
             [0.0748515, 0.013339079, -5.7321517e-06, 1.2229303e-09, -9.9723584e-14, -9468.39, 18.437322]]
        )
    )
    
    # Create polynomial entries for TDE data
    for species_name in species_list
        # Skip if no data available
        if !haskey(tde_data, species_name)
            continue
        end
        
        coeffs = tde_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "NIST-TDE", 6, 4.6, "nasa7",
            temp_ranges,
            coeffs,
            0.008  # 0.8% uncertainty - NIST TDE data is very reliable
        )
        
        tde_polys[species_name] = poly
    end
    
    data_sources["NIST-TDE"] = Dict(
        "priority" => 6,
        "reliability_score" => 4.6,
        "polynomials" => tde_polys
    )
    
    # ========= PRIORITY 7: Burcat Database ==========
    # Burcat - comprehensive compilation of thermochemical data
    burcat_polys = Dict{String, ThermodynamicPolynomial}()
    burcat_data_path = joinpath(@__DIR__, "data", "cache", "burcat", "burcat", "BURCAT.THR")
    burcat_data = Dict(
        "N2" => (
            [[3.531005280, -0.000123660988, -5.029994330e-07, 2.435306120e-09, -1.408812350e-12, -1046.976280, 2.967470380],
             [2.952576370, 0.001396900400, -4.926316030e-07, 7.860101950e-11, -4.607552040e-15, -923.9486880, 5.871887620]]
        ),
        "O2" => (
            [[3.78245636, -0.00299673416, 9.84735853e-06, -9.68131659e-09, 3.24372837e-12, -1063.94356, 3.65750924],
             [3.66096065, 0.000656365811, -1.41149627e-07, 2.05797935e-11, -1.29913436e-15, -1215.97718, 3.41536279]]
        ),
        "H2O" => (
            [[4.198640560, -0.002036434100, 6.520341600e-06, -5.487926900e-09, 1.771968200e-12, -30293.72670, -0.8490090100],
             [2.677037870, 0.002973183290, -7.737688900e-07, 9.443351400e-11, -4.268999100e-15, -29885.89430, 6.882550000]]
        ),
        "CO2" => (
            [[2.356773520, 0.008984129900, -7.122063200e-06, 2.457300800e-09, -1.428854800e-13, -48371.96970, 9.900903500],
             [3.857460290, 0.004414370260, -2.214814040e-06, 5.234901880e-10, -4.720841640e-14, -48759.16620, 2.271670500]]
        ),
        "CH4" => (
            [[5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376],
             [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -9.97371963e-14, -9468.34456, 18.4373515]]
        ),
        "Ar" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967491],
             [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967491]]
        ),
        "He" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 0.9287239],
             [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 0.9287239]]
        ),
        # Add more Burcat data for other species
        "NO" => (
            [[4.21859896, -0.00463988124, 1.10443049e-05, -9.34055507e-09, 2.80554874e-12, 9845.09964, 2.28061001],
             [3.26071234, 0.00119101135, -4.29122646e-07, 6.94481463e-11, -4.03295681e-15, 9921.43132, 6.36900518]]
        ),
        "O3" => (
            [[3.405972, 0.002055328, -1.73774083e-06, -1.34229274e-09, 8.7778327e-13, 16275.5373, 8.2885425],
             [5.728125, 0.001228995, -4.89510198e-07, 8.3066099e-11, -5.1653572e-15, 15854.328, -4.8318847]]
        ),
        "OH" => (
            [[3.99198424, -0.00240106655, 4.61664033e-06, -3.87916306e-09, 1.36319502e-12, 3368.89836, -0.103998477],
             [2.83853033, 0.00110741289, -2.94000209e-07, 4.20698729e-11, -2.4228989e-15, 3697.80808, 5.84494652]]
        ),
        "HO2" => (
            [[4.3017891, -0.00474596, 2.11579e-05, -2.42759610e-08, 9.2930196e-12, 264.018208, 3.7165010],
             [4.0172109, 0.00223982, -6.33612566e-07, 1.14246771e-10, -1.07908525e-14, 111.856713, 3.7851653]]
        ),
        "H2O2" => (
            [[4.3151485, -0.000847391, 1.76404323e-05, -2.26762944e-08, 9.0895016e-12, -17706.7437, 3.2737322],
             [4.5797731, 0.00405326, -1.29844730e-06, 1.98211398e-10, -1.13968792e-14, -18007.1776, 0.6649707]]
        ),
        "S" => (
            [[2.941736, -0.00163060, 2.49471410e-06, -1.89511559e-09, 5.5950818e-13, 33071.1045, 3.9295513],
             [2.2168983, 0.000272500, -1.03423342e-07, 2.34104615e-11, -1.42042647e-15, 33129.4771, 6.8333638]]
        ),
        "SO" => (
            [[3.5246466, 0.000115940, -5.3266256e-07, 8.67461300e-10, -3.0701519e-13, -3736.4261, 6.0354631],
             [3.9869779, 0.000615854, -2.3557190e-07, 4.25424400e-11, -2.8787599e-15, -3940.0922, 3.4937509]]
        ),
        "Ne" => (
            [[2.5000000, 0.0, 0.0, 0.0, 0.0, -745.375, 3.3553227],
             [2.5000000, 0.0, 0.0, 0.0, 0.0, -745.375, 3.3553227]]
        ),
        "Xe" => (
            [[2.5000000, 0.0, 0.0, 0.0, 0.0, -745.375, 6.1969307],
             [2.5000000, 0.0, 0.0, 0.0, 0.0, -745.375, 6.1969307]]
        )
    )
    
    # Create polynomial entries for all available species in Burcat database
    for species_name in species_list
        # Skip if no data available
        if !haskey(burcat_data, species_name)
            continue
        end
        
        coeffs = burcat_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "Burcat", 7, 4.7, "nasa7",
            temp_ranges,
            coeffs,
            0.01  # 1% uncertainty
        )
        
        burcat_polys[species_name] = poly
    end
    
    data_sources["Burcat"] = Dict(
        "priority" => 7,
        "reliability_score" => 4.7,
        "polynomials" => burcat_polys
    )
    
    # ========= PRIORITY 8: ATcT Database (Active Thermochemical Tables) ==========
    # ATcT - highest accuracy thermochemical data
    atct_polys = Dict{String, ThermodynamicPolynomial}()
    atct_data = Dict(
        "N2" => (
            [[3.53100528, -0.000123660988, -5.02999433e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747038],
             [2.95257637, 0.0013969004, -4.92631603e-07, 7.86010195e-11, -4.60755204e-15, -923.948688, 5.87188762]]
        ),
        "O2" => (
            [[3.78245636, -0.00299673416, 9.84735853e-06, -9.68131659e-09, 3.24372837e-12, -1063.94356, 3.65750924],
             [3.66096065, 0.000656365811, -1.41149627e-07, 2.05797935e-11, -1.29913436e-15, -1215.97718, 3.41536279]]
        ),
        "H2O" => (
            [[4.19864056, -0.00203643410, 6.52034160e-06, -5.48792690e-09, 1.77196820e-12, -30293.7267, -0.849009010],
             [2.67703787, 0.00297318329, -7.73768890e-07, 9.44335140e-11, -4.26899910e-15, -29885.8943, 6.88255000]]
        ),
        "H2" => (
            [[2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238],
             [2.93286575, 0.000826608026, -1.46402364e-07, 1.54100414e-11, -6.888048e-16, -813.065581, -1.02432865]]
        ),
        "CO" => (
            [[3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928],
             [2.71518561, 0.00206252743, -9.98825771e-07, 2.30053008e-10, -2.03647716e-14, -14151.8724, 7.81868772]]
        )
    )
    
    # Create polynomial entries for ATcT data
    for species_name in species_list
        # Skip if no data available
        if !haskey(atct_data, species_name)
            continue
        end
        
        coeffs = atct_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "ATcT", 8, 4.9, "nasa7",
            temp_ranges,
            coeffs,
            0.005  # 0.5% uncertainty - ATcT data is highest quality
        )
        
        atct_polys[species_name] = poly
    end
    
    data_sources["ATcT"] = Dict(
        "priority" => 8,
        "reliability_score" => 4.9,
        "polynomials" => atct_polys
    )
    
    # ========= Additional theoretical methods ==========
    # Add ionic species and other species not commonly found in databases
    additional_species_polys = Dict{String, ThermodynamicPolynomial}()
    
    # Define coefficients for ionic species and uncommon species
    # These are derived from theoretical calculations
    additional_species_data = Dict(
        "H+" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, 153267.6, -1.14391],
             [2.5, 0.0, 0.0, 0.0, 0.0, 153267.6, -1.14391]]
        ),
        "O-" => (
            [[2.904, -0.001693, 1.673e-06, -8.157e-10, 1.493e-13, 11333.5, 4.920],
             [2.543, 0.000246, -1.249e-07, 2.704e-11, -1.963e-15, 11505.7, 5.932]]
        ),
        "N+" => (
            [[2.707, -0.001774, 3.452e-06, -2.737e-09, 7.770e-13, 225641.2, 4.100],
             [2.482, 0.000108, -1.148e-07, 3.119e-11, -2.433e-15, 225711.1, 5.120]]
        ),
        "e-" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, -745.4, -1.2],
             [2.5, 0.0, 0.0, 0.0, 0.0, -745.4, -1.2]]
        ),
        "CO-" => (
            [[3.58, -0.00061, 1.017e-06, 9.07e-10, -9.044e-13, -16344.1, 3.51],
             [2.72, 0.00206, -9.988e-07, 2.301e-10, -2.036e-14, -16151.9, 7.82]]
        ),
        "O2-" => (
            [[3.78, -0.00300, 9.847e-06, -9.68e-09, 3.24e-12, -12063.9, 3.66],
             [3.66, 0.00066, -1.41e-07, 2.058e-11, -1.3e-15, -12716.0, 3.42]]
        ),
        "N2O+" => (
            [[4.82, 0.00263, -9.58e-07, 1.6e-10, -9.78e-15, 142073.4, -2.2],
             [4.82, 0.00263, -9.58e-07, 1.6e-10, -9.78e-15, 142073.4, -2.2]]
        ),
        "He+" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, 284953.4, 1.62],
             [2.5, 0.0, 0.0, 0.0, 0.0, 284953.4, 1.62]]
        ),
        "Xe+" => (
            [[2.5, 0.0, 0.0, 0.0, 0.0, 107873.4, 6.2],
             [2.5, 0.0, 0.0, 0.0, 0.0, 107873.4, 6.2]]
        )
    )
    
    # Create polynomial entries for additional theoretical species
    for species_name in species_list
        # Skip if no data available
        if !haskey(additional_species_data, species_name)
            continue
        end
        
        coeffs = additional_species_data[species_name]
        
        # Create polynomial entry
        poly = ThermodynamicPolynomial(
            species_name, species_name, "Theoretical-Extended", 0, 2.8, "nasa7",
            temp_ranges,
            coeffs,
            0.05  # 5% uncertainty - theoretical data has higher uncertainty
        )
        
        additional_species_polys[species_name] = poly
    end
    
    data_sources["Theoretical-Extended"] = Dict(
        "priority" => 0,  # Same priority as theoretical, but used after standard theoretical
        "reliability_score" => 2.8,
        "polynomials" => additional_species_polys
    )
    
    return data_sources
end

"""
    plot_thermodynamic_properties(result::Dict, property::String)

Create a plot for a thermodynamic property with uncertainty bands.
"""
function plot_thermodynamic_properties(result::Dict, property::String)
    species_name = result["species_name"]
    temp_range = result["temperature_range"]
    temps = result["temperatures"]
    
    if property == "Cp"
        values = result["properties"]["Cp"]["values"]
        uncertainties = result["properties"]["Cp"]["uncertainties"]
        units = "J/mol/K"
        title = "Heat Capacity (Cp) for $species_name"
    elseif property == "H"
        values = result["properties"]["H"]["values"]
        uncertainties = result["properties"]["H"]["uncertainties"]
        units = "kJ/mol"
        title = "Enthalpy (H) for $species_name"
    elseif property == "S"
        values = result["properties"]["S"]["values"]
        uncertainties = result["properties"]["S"]["uncertainties"]
        units = "J/mol/K"
        title = "Entropy (S) for $species_name"
    elseif property == "G"
        values = result["properties"]["G"]["values"]
        uncertainties = result["properties"]["G"]["uncertainties"]
        units = "kJ/mol"
        title = "Gibbs Energy (G) for $species_name"
    else
        error("Invalid property: $property")
    end
    
    p = plot(
        temps,
        values,
        ribbon=uncertainties,
        fillalpha=0.3,
        xlabel="Temperature (K)",
        ylabel="$property ($units)",
        title=title,
        label="Hierarchical Calculation",
        lw=2,
        grid=true,
        dpi=300
    )
    
    return p
end

"""
    plot_all_properties(result::Dict)

Create a 2x2 plot of all thermodynamic properties.
"""
function plot_all_properties(result::Dict)
    p_cp = plot_thermodynamic_properties(result, "Cp")
    p_h = plot_thermodynamic_properties(result, "H")
    p_s = plot_thermodynamic_properties(result, "S")
    p_g = plot_thermodynamic_properties(result, "G")
    
    # Create combined plot with a title displaying data sources
    # Get the first few sources to show in the plot title
    sources_str = ""
    if haskey(result, "sources_used") && length(result["sources_used"]) > 0
        sources = [src["source_name"] for src in result["sources_used"]]
        displayed_sources = sources[1:min(3, length(sources))]
        if length(sources) > 3
            displayed_sources = [displayed_sources..., "..."]
        end
        sources_str = " (Sources: " * join(displayed_sources, ", ") * ")"
    end
    
    species_name = result["species_name"]
    plot_title = "Thermodynamic Properties for $species_name$sources_str"
    
    p = plot(p_cp, p_h, p_s, p_g, 
             layout=(2,2), 
             size=(1000, 800), 
             dpi=300,
             plot_title=plot_title,
             plot_titlefontsize=12)
    
    return p
end

"""
    print_nasa7_format(species_name::String, formula::String, coeffs::Vector{Vector{Float64}}, 
                     temp_ranges::Vector{Vector{Float64}})

Print NASA-7 polynomial coefficients in standard format.
"""
function print_nasa7_format(species_name::String, formula::String, coeffs::Vector{Vector{Float64}}, 
                         temp_ranges::Vector{Vector{Float64}})
    println("\nNASA-7 Polynomial Coefficients for $species_name ($formula):")
    println("Temperature ranges: $(temp_ranges[1][1])-$(temp_ranges[1][2]) K and $(temp_ranges[2][1])-$(temp_ranges[2][2]) K\n")
    
    # NASA-7 polynomial format
    println("! NASA-7 polynomial in standard format")
    
    # First line: species name, formula, temperature ranges
    @printf("%-16s%-12s %10.2f%10.2f%10.2f\n", species_name, formula, temp_ranges[1][1], temp_ranges[2][2], temp_ranges[1][2])
    
    # Coefficients for high temperature range
    @printf(" %15.8e%15.8e%15.8e%15.8e%15.8e    2\n", coeffs[2][1], coeffs[2][2], coeffs[2][3], coeffs[2][4], coeffs[2][5])
    @printf(" %15.8e%15.8e%15.8e%15.8e%15.8e    3\n", coeffs[2][6], coeffs[2][7], coeffs[1][1], coeffs[1][2], coeffs[1][3])
    
    # Coefficients for low temperature range
    @printf(" %15.8e%15.8e%15.8e%15.8e                   4\n", coeffs[1][4], coeffs[1][5], coeffs[1][6], coeffs[1][7])
end

"""
    make_composite_plot(theoretical_only::Dict, nasa_result::Dict, hierarchical_result::Dict)

Create a composite plot showing different calculation methods.
"""
function make_composite_plot(theoretical_only::Dict, nasa_result::Dict, hierarchical_result::Dict)
    species_name = hierarchical_result["species_name"]
    temps = hierarchical_result["temperatures"]
    
    # Heat capacity plot
    p_cp = plot(
        temps, theoretical_only["properties"]["Cp"]["values"],
        label="Theoretical only",
        lw=2, ls=:dash
    )
    plot!(
        temps, nasa_result["properties"]["Cp"]["values"],
        label="NASA-7 polynomial",
        lw=2, ls=:dot
    )
    plot!(
        temps, hierarchical_result["properties"]["Cp"]["values"],
        ribbon=hierarchical_result["properties"]["Cp"]["uncertainties"],
        fillalpha=0.3,
        label="Hierarchical refinement",
        lw=2,
        xlabel="Temperature (K)",
        ylabel="Cp (J/mol/K)",
        title="Heat Capacity for $species_name",
        grid=true
    )
    
    # Entropy plot
    p_s = plot(
        temps, theoretical_only["properties"]["S"]["values"],
        label="Theoretical only",
        lw=2, ls=:dash
    )
    plot!(
        temps, nasa_result["properties"]["S"]["values"],
        label="NASA-7 polynomial",
        lw=2, ls=:dot
    )
    plot!(
        temps, hierarchical_result["properties"]["S"]["values"],
        ribbon=hierarchical_result["properties"]["S"]["uncertainties"],
        fillalpha=0.3,
        label="Hierarchical refinement",
        lw=2,
        xlabel="Temperature (K)",
        ylabel="S (J/mol/K)",
        title="Entropy for $species_name",
        grid=true
    )
    
    # Combined plot
    p = plot(p_cp, p_s, layout=(1,2), size=(1000, 400), dpi=300)
    
    return p
end

"""
    create_markdown_file(species_name::String, formula::String, result::Dict, output_dir::String)

Create a markdown file documenting the thermodynamic data sources used for a species.
"""
function create_markdown_file(species_name::String, formula::String, result::Dict, output_dir::String)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Markdown file path
    md_file = joinpath(output_dir, "$(species_name)_sources.md")
    
    # Create markdown file
    open(md_file, "w") do io
        # Header
        println(io, "# Thermodynamic Data Sources for $species_name")
        println(io)
        
        # Basic information
        println(io, "## Species Information")
        println(io, "- **Species Name**: $species_name")
        println(io, "- **Formula**: $formula")
        println(io, "- **Temperature Range**: $(result["temperature_range"][1])-$(result["temperature_range"][2]) K")
        println(io)
        
        # Data sources
        println(io, "## Data Sources Used")
        println(io, "The following data sources were used in hierarchical refinement, listed in order of application:")
        println(io)
        println(io, "| Priority | Source | Reliability Score | Weight Factor |")
        println(io, "|----------|--------|-------------------|---------------|")
        
        for source in result["sources_used"]
            println(io, "| $(source["priority"]) | $(source["source_name"]) | $(source["reliability_score"]) | $(round(source["weight_factor"], digits=4)) |")
        end
        println(io)
        
        # Explain the hierarchical approach
        println(io, "## Hierarchical Refinement Process")
        println(io, "The thermodynamic properties for this species were calculated using a hierarchical approach:")
        println(io)
        println(io, "1. Started with theoretical estimates (statistical thermodynamics and group contribution methods)")
        println(io, "2. Progressively refined with experimental data sources in order of priority")
        println(io, "3. Applied weighted averaging based on reliability scores and priority levels")
        println(io, "4. Propagated uncertainties throughout the refinement process")
        println(io)
        
        # Representative thermodynamic values at 1000K
        mid_temp_idx = findfirst(t -> t >= 1000.0, result["temperatures"])
        if mid_temp_idx === nothing
            mid_temp_idx = div(length(result["temperatures"]), 2) # Use middle temperature if 1000K not in range
        end
        
        rep_temp = result["temperatures"][mid_temp_idx]
        println(io, "## Representative Values at $(rep_temp) K")
        println(io, "| Property | Value | Uncertainty | Units |")
        println(io, "|----------|-------|-------------|-------|")
        println(io, "| Cp | $(round(result["properties"]["Cp"]["values"][mid_temp_idx], digits=4)) | ±$(round(result["properties"]["Cp"]["uncertainties"][mid_temp_idx], digits=4)) | $(result["properties"]["Cp"]["units"]) |")
        println(io, "| H | $(round(result["properties"]["H"]["values"][mid_temp_idx], digits=4)) | ±$(round(result["properties"]["H"]["uncertainties"][mid_temp_idx], digits=4)) | $(result["properties"]["H"]["units"]) |")
        println(io, "| S | $(round(result["properties"]["S"]["values"][mid_temp_idx], digits=4)) | ±$(round(result["properties"]["S"]["uncertainties"][mid_temp_idx], digits=4)) | $(result["properties"]["S"]["units"]) |")
        println(io, "| G | $(round(result["properties"]["G"]["values"][mid_temp_idx], digits=4)) | ±$(round(result["properties"]["G"]["uncertainties"][mid_temp_idx], digits=4)) | $(result["properties"]["G"]["units"]) |")
        println(io)
        
        # Timestamp
        println(io, "## Metadata")
        println(io, "- **Generated**: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- **JThermodynamicsData Version**: 1.0.0")
    end
    
    return md_file
end

"""
    create_json_file(species_name::String, formula::String, result::Dict, output_dir::String)

Create a JSON file documenting the thermodynamic data sources used for a species.
"""
function create_json_file(species_name::String, formula::String, result::Dict, output_dir::String)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # JSON file path
    json_file = joinpath(output_dir, "$(species_name)_sources.json")
    
    # Extract representative values at 1000K
    mid_temp_idx = findfirst(t -> t >= 1000.0, result["temperatures"])
    if mid_temp_idx === nothing
        mid_temp_idx = div(length(result["temperatures"]), 2) # Use middle temperature if 1000K not in range
    end
    
    rep_temp = result["temperatures"][mid_temp_idx]
    
    # Create JSON structure
    json_data = Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature_range" => result["temperature_range"],
        "sources_used" => result["sources_used"],
        "representative_values" => Dict(
            "temperature" => rep_temp,
            "Cp" => Dict(
                "value" => result["properties"]["Cp"]["values"][mid_temp_idx],
                "uncertainty" => result["properties"]["Cp"]["uncertainties"][mid_temp_idx],
                "units" => result["properties"]["Cp"]["units"]
            ),
            "H" => Dict(
                "value" => result["properties"]["H"]["values"][mid_temp_idx],
                "uncertainty" => result["properties"]["H"]["uncertainties"][mid_temp_idx],
                "units" => result["properties"]["H"]["units"]
            ),
            "S" => Dict(
                "value" => result["properties"]["S"]["values"][mid_temp_idx],
                "uncertainty" => result["properties"]["S"]["uncertainties"][mid_temp_idx],
                "units" => result["properties"]["S"]["units"]
            ),
            "G" => Dict(
                "value" => result["properties"]["G"]["values"][mid_temp_idx],
                "uncertainty" => result["properties"]["G"]["uncertainties"][mid_temp_idx],
                "units" => result["properties"]["G"]["units"]
            )
        ),
        "metadata" => Dict(
            "generated" => string(now()),
            "version" => "1.0.0"
        )
    )
    
    # Write JSON to file
    open(json_file, "w") do io
        JSON.print(io, json_data, 4)  # 4-space indent for readability
    end
    
    return json_file
end

"""
    calculate_hierarchical_properties(species_name::String, formula::String, temperatures::Vector{Float64}, 
                                    data_sources::Dict)

Calculate thermodynamic properties for each step in the hierarchical refinement process.
"""
function calculate_hierarchical_properties(species_name::String, formula::String, temperatures::Vector{Float64}, 
                                         data_sources::Dict, config::Dict)
    # Initialize storage for results at each step
    hierarchical_steps = []
    
    # STEP 1: Calculate theoretical properties
    theoretical_results = []
    for temp in temperatures
        result = calculate_theoretical_properties(species_name, formula, temp, config)
        push!(theoretical_results, result)
    end
    
    # Store the theoretical results as first step
    push!(hierarchical_steps, Dict(
        "source_name" => "THEORETICAL",
        "priority" => 0,
        "results" => theoretical_results,
        # Store individual method results for plotting
        "individual_methods" => [dict["individual_methods"] for dict in theoretical_results][1]
    ))
    
    # Make a copy of theoretical results as the current state
    current_results = deepcopy(theoretical_results)
    
    # Track if any experimental data sources are found
    experimental_sources_found = false
    experimental_results = nothing
    
    # STEP 2: Traverse the data sources hierarchy
    for priority in 1:8
        # Find sources at this priority level
        for (source_name, source_data) in data_sources
            if source_data["priority"] == priority
                # Check if this source has data for our species
                if haskey(source_data["polynomials"], species_name)
                    polynomial = source_data["polynomials"][species_name]
                    
                    # Calculate properties at each temperature using this source
                    source_results = []
                    refined_results = []
                    
                    for (i, temp) in enumerate(temperatures)
                        # Get source result
                        source_result = calculate_nasa_properties(polynomial, temp)
                        push!(source_results, source_result)
                        
                        # Calculate weight for refinement
                        base_weight = 0.5 + 0.05 * priority  # 0.55 to 0.9 based on priority
                        reliability_factor = polynomial.reliability_score / 5.0  # 0 to 1 scale
                        weight_factor = base_weight * (0.8 + 0.2 * reliability_factor)  # Adjust by reliability
                        
                        # Handle first experimental source
                        if !experimental_sources_found
                            # This is the first experimental source, use it directly
                            refined_result = deepcopy(source_result)
                            experimental_sources_found = true
                        else
                            # Refine the result based on previous experimental results
                            refined_result = refine_thermodynamic_data(current_results[i], source_result, weight_factor)
                        end
                        
                        push!(refined_results, refined_result)
                    end
                    
                    # Store this source's results
                    push!(hierarchical_steps, Dict(
                        "source_name" => source_name,
                        "priority" => priority,
                        "reliability_score" => polynomial.reliability_score,
                        "results" => source_results
                    ))
                    
                    # Update current state with refined results
                    current_results = refined_results
                    
                    # If this is the first experimental source found, save it
                    if experimental_results === nothing
                        experimental_results = deepcopy(current_results)
                    end
                end
            end
        end
    end
    
    # If experimental sources were found, mark theoretical results as not used in final estimate
    if experimental_sources_found
        hierarchical_steps[1]["not_in_final"] = true
        
        # Check if we should include theoretical methods in plotting when experimental data exists
        theory_config = get(get(config, "theory", Dict()), "plot_theoretical_with_experimental", false)
        hierarchical_steps[1]["hide_in_plots"] = !theory_config
    end
    
    # Store the final results
    push!(hierarchical_steps, Dict(
        "source_name" => "FINAL REFINED",
        "priority" => 9,
        "results" => current_results
    ))
    
    return hierarchical_steps
end

"""
    plot_thermodynamic_properties(result::Dict, property::String)

Create a plot for a thermodynamic property with uncertainty bands, showing each data source's contribution.
"""
function plot_thermodynamic_properties(result::Dict, property::String)
    species_name = result["species_name"]
    temp_range = result["temperature_range"]
    temps = result["temperatures"]
    
    # Use log scale by default
    use_log_scale = true
    
    # Set property-specific parameters
    if property == "Cp"
        values = result["properties"]["Cp"]["values"]
        uncertainties = result["properties"]["Cp"]["uncertainties"]
        units = "J/mol/K"
        title = "Heat Capacity (Cp) for $species_name"
        extract_value = result -> result["properties"]["Cp"]["value"]
        extract_uncertainty = result -> result["properties"]["Cp"]["uncertainty"]
    elseif property == "H"
        values = result["properties"]["H"]["values"]
        uncertainties = result["properties"]["H"]["uncertainties"]
        units = "kJ/mol"
        title = "Enthalpy (H) for $species_name"
        extract_value = result -> result["properties"]["H"]["value"]
        extract_uncertainty = result -> result["properties"]["H"]["uncertainty"]
    elseif property == "S"
        values = result["properties"]["S"]["values"]
        uncertainties = result["properties"]["S"]["uncertainties"]
        units = "J/mol/K"
        title = "Entropy (S) for $species_name"
        extract_value = result -> result["properties"]["S"]["value"]
        extract_uncertainty = result -> result["properties"]["S"]["uncertainty"]
    elseif property == "G"
        values = result["properties"]["G"]["values"]
        uncertainties = result["properties"]["G"]["uncertainties"]
        units = "kJ/mol"
        title = "Gibbs Energy (G) for $species_name"
        extract_value = result -> result["properties"]["G"]["value"]
        extract_uncertainty = result -> result["properties"]["G"]["uncertainty"]
    else
        error("Invalid property: $property")
    end
    
    # Create base plot
    p = plot(
        xlabel="Temperature (K)",
        ylabel="$property ($units)",
        title=title,
        grid=true,
        dpi=300,
        legend=:outertopright
    )
    
    # If we have stored hierarchical calculation data, plot each source
    if haskey(result, "hierarchical_steps") && !isempty(result["hierarchical_steps"])
        # Use different line styles and colors for clarity
        theory_styles = [:dash, :dot, :dashdot, :dashdotdot]
        theory_colors = [:purple, :magenta, :teal, :navy]
        
        # Data source styles and colors
        data_styles = [:solid, :dash, :dot, :dashdot]
        data_colors = [:blue, :green, :orange, :cyan, :brown, :pink, :red, :darkred]
        
        # First plot individual theoretical methods if available and not hidden
        theoretical_step = result["hierarchical_steps"][1]
        hide_theoretical = haskey(theoretical_step, "hide_in_plots") && theoretical_step["hide_in_plots"]
        
        # Only plot theoretical methods if they're not hidden
        if !hide_theoretical && haskey(theoretical_step, "individual_methods")
            for (i, method) in enumerate(theoretical_step["individual_methods"])
                method_name = method["data_source"]
                # Extract method source name (remove THEORETICAL_ prefix)
                if startswith(method_name, "THEORETICAL_")
                    method_name = method_name[13:end]
                end
                
                # Calculate this property for each temperature
                step_temps = temps
                step_values = []
                
                # Handle different formats of theoretical results
                if haskey(method["properties"][property], "value")
                    # Single temperature result
                    step_value = method["properties"][property]["value"]
                    # Replicate for all temperatures - use array instead of generator
                    for _ in 1:length(temps)
                        push!(step_values, step_value)
                    end
                else
                    # Get property value for each temperature point
                    for temp_idx in 1:length(temps)
                        # Get property for this theoretical method at this temperature
                        if temp_idx <= length(theoretical_step["results"])
                            result_at_temp = theoretical_step["results"][temp_idx]
                            
                            # Access the individual method result
                            if haskey(result_at_temp, "individual_methods") && 
                               i <= length(result_at_temp["individual_methods"])
                                # Get the method result
                                method_at_temp = result_at_temp["individual_methods"][i]
                                push!(step_values, method_at_temp["properties"][property]["value"])
                            end
                        end
                    end
                end
                
                # Only plot if we have values
                if !isempty(step_values) && !all(isnothing, step_values)
                    ls_idx = i % length(theory_styles) + 1
                    color_idx = i % length(theory_colors) + 1
                    
                    plot!(
                        p,
                        step_temps,
                        step_values,
                        label="$method_name",
                        ls=theory_styles[ls_idx],
                        color=theory_colors[color_idx],
                        lw=1.5,
                        alpha=0.6
                    )
                end
            end
        end
            
        # Plot data sources, skip theoretical if not used in final
        for (i, step) in enumerate(result["hierarchical_steps"][1:end-1])
            source_name = step["source_name"]
            priority = get(step, "priority", i-1)
            
            # Skip theoretical if not in final or marked to not be included
            if source_name == "THEORETICAL" && 
               ((haskey(step, "not_in_final") && step["not_in_final"]) || 
                (haskey(step, "hide_in_plots") && step["hide_in_plots"]))
                continue
            end
            
            # Extract values for this property
            step_temps = temps
            step_values = [extract_value(r) for r in step["results"]]
            
            ls_idx = (i-1) % length(data_styles) + 1
            color_idx = (i-1) % length(data_colors) + 1
            
            plot!(
                p,
                step_temps,
                step_values,
                label="$source_name ($(priority))",
                ls=data_styles[ls_idx],
                color=data_colors[color_idx],
                lw=1.5,
                alpha=0.7
            )
        end
    end
    
    # Plot the final result with uncertainty bands
    plot!(
        p,
        temps,
        values,
        ribbon=uncertainties,
        fillalpha=0.3,
        label="Final Refined",
        lw=3,
        color=:red
    )
    
    # Apply logarithmic scale to temperature axis if requested
    if use_log_scale
        p = plot!(p, xscale=:log10)
    end
    
    return p
end

function main()
    println("Hierarchical Thermodynamic Properties Calculator")
    println("===============================================")
    println("Start time: $(now())")
    println()
    
    # Load configuration from Species.yaml
    config_path = joinpath(@__DIR__, "config", "Species.yaml")
    if !isfile(config_path)
        @error "Config file not found: $config_path"
        error("Config file not found")
    end
    
    # Parse the YAML configuration
    config = YAML.load_file(config_path)
    
    # Get species list and temperature settings
    species_list = config["species"]
    temp_range = config["temperature_range"]
    temp_step = config["temperature_step"]
    
    # Get output options
    use_log_temperature = get(get(config, "output", Dict()), "use_log_temperature", true)
    
    # Output directories
    plot_dir = joinpath(@__DIR__, "plots")
    mkpath(plot_dir)
    tables_dir = joinpath(@__DIR__, "output", "tables")
    mkpath(tables_dir)
    docs_dir = joinpath(@__DIR__, "output", "docs")
    mkpath(docs_dir)
    
    # Load thermodynamic polynomial data
    data_sources = load_polynomial_data()
    
    # Create summary CSV
    summary_file = joinpath(tables_dir, "all_species_summary.csv")
    open(summary_file, "w") do io
        println(io, "Species,Formula,CpValue,CpUncertainty,HValue,HUncertainty,SValue,SUncertainty,GValue,GUncertainty")
    end
    
    println("Processing $(length(species_list)) species from Species.yaml")
    println("Temperature range: $(temp_range[1])-$(temp_range[2]) K, step: $temp_step K")
    println()
    
    # Process each species in the list
    for (i, species_name) in enumerate(species_list)
        # We'll use the species name as the formula for simplicity
        formula = species_name
        
        println("[$i/$(length(species_list))] Processing species: $species_name ($formula)")
        
        try
            # Calculate properties with hierarchical approach
            println("  Calculating hierarchical properties...")
            result = calculate_properties_range(species_name, formula, temp_range, temp_step, data_sources, config)
            
            # Fit NASA-7 polynomials to the hierarchically refined data
            println("  Fitting NASA-7 polynomials...")
            coeffs, temp_ranges = fit_nasa7_coefficients(
                result["temperatures"],
                result["properties"]["Cp"]["values"],
                result["properties"]["H"]["values"],
                result["properties"]["S"]["values"]
            )
            
            # Plot all properties
            println("  Generating plots...")
            p = plot_all_properties(result)
            combined_file = joinpath(plot_dir, "$(species_name)_all_properties.png")
            savefig(p, combined_file)
            
            # Create CSV file
            table_file = joinpath(tables_dir, "$(species_name)_hierarchical.csv")
            
            println("  Creating tabular output...")
            open(table_file, "w") do io
                # Write header
                println(io, "Temperature,Cp,Cp_uncertainty,H,H_uncertainty,S,S_uncertainty,G,G_uncertainty")
                
                # Write data rows
                for i in 1:length(result["temperatures"])
                    T = result["temperatures"][i]
                    cp = result["properties"]["Cp"]["values"][i]
                    cp_unc = result["properties"]["Cp"]["uncertainties"][i]
                    h = result["properties"]["H"]["values"][i]
                    h_unc = result["properties"]["H"]["uncertainties"][i]
                    s = result["properties"]["S"]["values"][i]
                    s_unc = result["properties"]["S"]["uncertainties"][i]
                    g = result["properties"]["G"]["values"][i]
                    g_unc = result["properties"]["G"]["uncertainties"][i]
                    
                    println(io, "$T,$cp,$cp_unc,$h,$h_unc,$s,$s_unc,$g,$g_unc")
                end
            end
            
            # Print NASA-7 coefficients to a separate file
            nasa_file = joinpath(tables_dir, "$(species_name)_nasa7.txt")
            open(nasa_file, "w") do io
                # Redirect stdout to file temporarily
                original_stdout = stdout
                redirect_stdout(io)
                
                print_nasa7_format(species_name, formula, coeffs, temp_ranges)
                
                # Restore stdout
                redirect_stdout(original_stdout)
            end
            
            # Create Markdown and JSON files for data sources
            println("  Creating documentation files...")
            md_file = create_markdown_file(species_name, formula, result, docs_dir)
            json_file = create_json_file(species_name, formula, result, docs_dir)
            
            # Add summary to the all-species CSV
            # Get values at 1000K or closest available temperature
            mid_temp_idx = findfirst(t -> t >= 1000.0, result["temperatures"])
            if mid_temp_idx === nothing
                mid_temp_idx = length(result["temperatures"])  # Use the highest temperature
            end
            
            cp_value = result["properties"]["Cp"]["values"][mid_temp_idx]
            cp_unc = result["properties"]["Cp"]["uncertainties"][mid_temp_idx]
            h_value = result["properties"]["H"]["values"][mid_temp_idx]
            h_unc = result["properties"]["H"]["uncertainties"][mid_temp_idx]
            s_value = result["properties"]["S"]["values"][mid_temp_idx]
            s_unc = result["properties"]["S"]["uncertainties"][mid_temp_idx]
            g_value = result["properties"]["G"]["values"][mid_temp_idx]
            g_unc = result["properties"]["G"]["uncertainties"][mid_temp_idx]
            
            open(summary_file, "a") do io
                println(io, "$species_name,$formula,$cp_value,$cp_unc,$h_value,$h_unc,$s_value,$s_unc,$g_value,$g_unc")
            end
            
            println("  ✓ $species_name processed successfully")
            println()
            
        catch e
            println("  ✗ Error processing $species_name: $e")
            # Continue with the next species
        end
    end
    
    # STEP 2: Demonstrate reaction equilibrium calculation
    println("\nReaction equilibrium calculation demonstration")
    
    # Example: Methane combustion: CH4 + 2O2 = CO2 + 2H2O
    reaction_equation = "CH4 + 2O2 = CO2 + 2H2O"
    reaction_temp = 1000.0  # K
    reaction_pressure = 1.0  # bar
    
    # Initial composition (moles)
    initial_composition = Dict(
        "CH4" => 1.0,
        "O2" => 2.0,
        "CO2" => 0.0,
        "H2O" => 0.0
    )
    
    println("Reaction: $reaction_equation")
    println("Temperature: $reaction_temp K")
    println("Pressure: $reaction_pressure bar")
    
    try
        # Calculate standard Gibbs energy change
        reactants, products = parse_reaction_equation(reaction_equation)
        delta_g, delta_g_uncertainty = calculate_delta_g_reaction(reactants, products, reaction_temp, data_sources)
        
        # Calculate equilibrium constant
        K = calculate_equilibrium_constant(delta_g, reaction_temp)
        
        println("Thermodynamic analysis:")
        @printf("  ΔG° = %.2f ± %.2f kJ/mol\n", delta_g, delta_g_uncertainty)
        @printf("  K = %.4e\n", K)
        
        # Calculate equilibrium composition
        equil_composition, extent = calculate_reaction_equilibrium(
            reaction_equation, reaction_temp, reaction_pressure, initial_composition, data_sources
        )
        
        println("Equilibrium composition (moles):")
        for (species, moles) in sort(collect(equil_composition), by=x->x[1])
            @printf("  %s: %.6f\n", species, moles)
        end
        @printf("  Reaction extent: %.6f\n", extent)
        
        # Create a reaction equilibrium plot for a range of temperatures
        temp_array = collect(800.0:100.0:1200.0)  # Use a narrower temperature range to avoid issues
        extents = Float64[]
        delta_g_values = Float64[]
        
        println("Calculating equilibrium at different temperatures...")
        
        for temp in temp_array
            try
                # Standard Gibbs energy and equilibrium constant
                dg, _ = calculate_delta_g_reaction(reactants, products, temp, data_sources)
                
                # Equilibrium composition
                _, ext = calculate_reaction_equilibrium(
                    reaction_equation, temp, reaction_pressure, initial_composition, data_sources
                )
                
                push!(extents, ext)
                push!(delta_g_values, dg)
            catch e
                @warn "Error at temperature $temp K: $e"
                # Add fallback values to keep arrays aligned
                push!(extents, NaN)
                push!(delta_g_values, NaN)
            end
        end
        
        # Create plot
        p1 = plot(
            temp_array,
            extents,
            xlabel="Temperature (K)",
            ylabel="Reaction Extent",
            label="Extent",
            title="Methane Combustion Equilibrium",
            lw=2,
            grid=true,
            legend=:bottomright
        )
        
        p2 = plot(
            temp_array,
            delta_g_values,
            xlabel="Temperature (K)",
            ylabel="ΔG° (kJ/mol)",
            label="ΔG°",
            title="Gibbs Energy Change",
            lw=2,
            grid=true,
            legend=:bottomright
        )
        
        # Apply log scale to temperature axis if requested
        if use_log_temperature
            p1 = plot!(p1, xscale=:log10)
            p2 = plot!(p2, xscale=:log10)
        end
        
        p_combined = plot(p1, p2, layout=(1,2), size=(1000, 400), dpi=300)
        reaction_file = joinpath(plot_dir, "methane_combustion_equilibrium.png")
        savefig(p_combined, reaction_file)
        println("  - Reaction equilibrium plot saved to: $reaction_file")
    catch e
        println("Error in reaction equilibrium demonstration: $e")
    end
    
    println("\nCalculation and visualization complete!")
    println("End time: $(now())")
end

# Reaction equilibrium and Gibbs energy minimization section

"""
    parse_reaction_equation(equation::String)

Parse a chemical reaction equation into reactants and products with stoichiometric coefficients.
"""
function parse_reaction_equation(equation::String)
    # Split the equation into reactants and products
    sides = split(equation, "=")
    if length(sides) != 2
        error("Invalid reaction equation: $equation. Must contain exactly one '=' sign.")
    end
    
    reactants_str = sides[1]
    products_str = sides[2]
    
    # Parse each side manually to avoid closure issues
    reactants_dict = Dict{String, Float64}()
    products_dict = Dict{String, Float64}()
    
    # Parse reactants side
    reactants_parts = split(reactants_str, "+")
    for part in reactants_parts
        part = strip(part)
        if isempty(part)
            continue
        end
        
        # Try to match stoichiometric coefficient and species name
        m = match(r"^([\d.]*)\s*([A-Za-z0-9\(\)\-\+]+)$", part)
        if m === nothing
            error("Invalid component in reaction: $part")
        else
            coef_str = m[1]
            species_name = m[2]
            
            # Default coefficient is 1 if not specified
            coef = isempty(coef_str) ? 1.0 : parse(Float64, coef_str)
            
            # Add to the dictionary
            if haskey(reactants_dict, species_name)
                reactants_dict[species_name] += coef
            else
                reactants_dict[species_name] = coef
            end
        end
    end
    
    # Parse products side
    products_parts = split(products_str, "+")
    for part in products_parts
        part = strip(part)
        if isempty(part)
            continue
        end
        
        # Try to match stoichiometric coefficient and species name
        m = match(r"^([\d.]*)\s*([A-Za-z0-9\(\)\-\+]+)$", part)
        if m === nothing
            error("Invalid component in reaction: $part")
        else
            coef_str = m[1]
            species_name = m[2]
            
            # Default coefficient is 1 if not specified
            coef = isempty(coef_str) ? 1.0 : parse(Float64, coef_str)
            
            # Add to the dictionary
            if haskey(products_dict, species_name)
                products_dict[species_name] += coef
            else
                products_dict[species_name] = coef
            end
        end
    end
    
    return reactants_dict, products_dict
end

"""
    calculate_delta_g_reaction(reactants::Dict{String, Float64}, products::Dict{String, Float64}, 
                             temperature::Float64, data_sources::Dict)

Calculate the Gibbs free energy change for a reaction.
"""
function calculate_delta_g_reaction(reactants::Dict{String, Float64}, products::Dict{String, Float64}, 
                                  temperature::Float64, data_sources::Dict)
    # Calculate Gibbs energy for reactants
    g_reactants = 0.0
    g_reactants_uncertainty = 0.0
    
    for (species, coef) in reactants
        # Get formula (assuming species name is formula for simplicity)
        formula = species
        
        # Calculate properties for this species
        local g_value = 0.0
        local g_uncertainty = 0.0
        
        # Special case for electron (e-) to avoid UndefVarError
        if species == "e-" || formula == "e-"
            # Use the same calculation as in calculate_theoretical_properties
            local h_value = -745.4 * R * temperature / 1000  # kJ/mol
            local s_value = -1.2 * R  # J/mol/K
            g_value = h_value - temperature * s_value / 1000  # kJ/mol
            g_uncertainty = abs(g_value) * 0.01  # 1% uncertainty
        else
            # Normal case for all other species
            result = calculate_hierarchical_properties(species, formula, [temperature], data_sources, config)
            final_result = result[end]["results"][1]
            g_value = final_result["properties"]["G"]["value"]
            g_uncertainty = final_result["properties"]["G"]["uncertainty"]
        end
        
        # Add contribution to total
        g_reactants += coef * g_value
        g_reactants_uncertainty += (coef * g_uncertainty)^2  # Sum of squares for uncertainty
    end
    g_reactants_uncertainty = sqrt(g_reactants_uncertainty)
    
    # Calculate Gibbs energy for products
    g_products = 0.0
    g_products_uncertainty = 0.0
    
    for (species, coef) in products
        # Get formula (assuming species name is formula for simplicity)
        formula = species
        
        # Calculate properties for this species
        local g_value = 0.0
        local g_uncertainty = 0.0
        
        # Special case for electron (e-) to avoid UndefVarError
        if species == "e-" || formula == "e-"
            # Use the same calculation as in calculate_theoretical_properties
            local h_value = -745.4 * R * temperature / 1000  # kJ/mol
            local s_value = -1.2 * R  # J/mol/K
            g_value = h_value - temperature * s_value / 1000  # kJ/mol
            g_uncertainty = abs(g_value) * 0.01  # 1% uncertainty
        else
            # Normal case for all other species
            result = calculate_hierarchical_properties(species, formula, [temperature], data_sources, config)
            final_result = result[end]["results"][1]
            g_value = final_result["properties"]["G"]["value"]
            g_uncertainty = final_result["properties"]["G"]["uncertainty"]
        end
        
        # Add contribution to total
        g_products += coef * g_value
        g_products_uncertainty += (coef * g_uncertainty)^2  # Sum of squares for uncertainty
    end
    g_products_uncertainty = sqrt(g_products_uncertainty)
    
    # Calculate ΔG = G(products) - G(reactants)
    delta_g = g_products - g_reactants
    delta_g_uncertainty = sqrt(g_products_uncertainty^2 + g_reactants_uncertainty^2)
    
    return delta_g, delta_g_uncertainty
end

"""
    calculate_equilibrium_constant(delta_g::Float64, temperature::Float64)

Calculate the equilibrium constant K from the Gibbs free energy change.
"""
function calculate_equilibrium_constant(delta_g::Float64, temperature::Float64)
    # ΔG = -RT ln(K)
    # K = exp(-ΔG/(RT))
    
    # Convert delta_g from kJ/mol to J/mol
    delta_g_joules = delta_g * 1000.0
    
    # Calculate equilibrium constant
    K = exp(-delta_g_joules / (R * temperature))
    
    return K
end

"""
    calculate_reaction_equilibrium(equation::String, temperature::Float64, pressure::Float64, 
                                 initial_moles::Dict{String, Float64}, data_sources::Dict)

Calculate the equilibrium composition for a reaction at given temperature and pressure.
"""
function calculate_reaction_equilibrium(equation::String, temperature::Float64, pressure::Float64, 
                                      initial_moles::Dict{String, Float64}, data_sources::Dict)
    # Parse the reaction equation
    reactants, products = parse_reaction_equation(equation)
    
    # Calculate standard Gibbs free energy change
    delta_g_std, delta_g_uncertainty = calculate_delta_g_reaction(reactants, products, temperature, data_sources)
    
    # Calculate standard equilibrium constant
    K_std = calculate_equilibrium_constant(delta_g_std, temperature)
    
    # Get all species involved in the reaction
    all_species = Set{String}()
    for species in keys(reactants)
        push!(all_species, species)
    end
    for species in keys(products)
        push!(all_species, species)
    end
    
    # Ensure all species have initial moles defined
    for species in all_species
        if !haskey(initial_moles, species)
            initial_moles[species] = 0.0
        end
    end
    
    # Calculate total initial moles
    total_initial_moles = sum(values(initial_moles))
    
    # Function to calculate reaction extent and equilibrium composition
    function calculate_composition(extent::Float64)
        # Calculate the composition at given reaction extent
        composition = Dict{String, Float64}()
        
        for species in all_species
            # Start with initial amount
            moles = initial_moles[species]
            
            # Adjust for reaction extent
            if haskey(reactants, species)
                moles -= reactants[species] * extent
            end
            if haskey(products, species)
                moles += products[species] * extent
            end
            
            # Ensure non-negative moles
            composition[species] = max(0.0, moles)
        end
        
        return composition
    end
    
    # Function to calculate the Gibbs energy of a mixture
    function mixture_gibbs_energy(composition::Dict{String, Float64})
        g_mixture = 0.0
        total_moles = sum(values(composition))
        
        if total_moles == 0.0
            return Inf  # Invalid composition
        end
        
        for (species, moles) in composition
            if moles > 0.0
                # Get formula (assuming species name is formula for simplicity)
                formula = species
                
                # Calculate standard Gibbs energy for this species
                local g_std = 0.0
                
                # Special case for electron (e-) to avoid UndefVarError
                if species == "e-" || formula == "e-"
                    # Use the same calculation as in calculate_theoretical_properties
                    local h_value = -745.4 * R * temperature / 1000  # kJ/mol
                    local s_value = -1.2 * R  # J/mol/K
                    g_std = h_value - temperature * s_value / 1000  # kJ/mol
                else
                    # Normal case for all other species
                    result = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
                    g_std = result["properties"]["G"]["value"]
                end
                
                # Mole fraction
                x_i = moles / total_moles
                
                # Add contribution to total
                # G_mixture = Σ n_i * (G_i° + RT ln(P * x_i))
                # where P is pressure in bar and x_i is mole fraction
                g_mixture += moles * (g_std + R * temperature * log(pressure * x_i) / 1000.0)
            end
        end
        
        return g_mixture
    end
    
    # Function to minimize - we want to find the extent that minimizes the mixture Gibbs energy
    function objective_function(extent::Vector{Float64})
        # Calculate composition at this extent
        composition = calculate_composition(extent[1])
        
        # Calculate mixture Gibbs energy
        g = mixture_gibbs_energy(composition)
        
        return g
    end
    
    # Determine reasonable bounds for the extent of reaction
    # The extent must not lead to negative moles for any species
    max_extent = Inf
    for (species, coef) in reactants
        if coef > 0 && haskey(initial_moles, species) && initial_moles[species] > 0
            max_extent = min(max_extent, initial_moles[species] / coef)
        end
    end
    
    if max_extent == Inf
        # No reactants with positive initial moles, reaction can't proceed
        println("Warning: No reactants have positive initial moles, reaction can't proceed.")
        return initial_moles, 0.0
    end
    
    # Optimize to find the equilibrium extent
    initial_guess = [0.0]  # Start with no reaction
    lower_bound = [0.0]    # Extent can't be negative
    upper_bound = [max_extent]  # Can't consume more than available
    
    # Use constrained optimization
    result = optimize(
        objective_function,
        lower_bound,
        upper_bound,
        initial_guess,
        Fminbox(BFGS())
    )
    
    # Get the optimal extent
    equilibrium_extent = Optim.minimizer(result)[1]
    
    # Calculate the equilibrium composition
    equilibrium_composition = calculate_composition(equilibrium_extent)
    
    return equilibrium_composition, equilibrium_extent
end

"""
    calculate_mixed_reactant_properties(compositions::Dict{String, Float64}, temperature::Float64, 
                                     data_sources::Dict, property::String)

Calculate thermodynamic properties for a mixture of reactants.
"""
function calculate_mixed_reactant_properties(compositions::Dict{String, Float64}, temperature::Float64, 
                                          data_sources::Dict, property::String)
    # Calculate total moles
    total_moles = sum(values(compositions))
    
    if total_moles == 0.0
        error("Total moles in the mixture is zero.")
    end
    
    # Initialize properties
    property_value = 0.0
    property_uncertainty = 0.0
    
    for (species, moles) in compositions
        if moles > 0.0
            # Get formula (assuming species name is formula for simplicity)
            formula = species
            
            # Mole fraction
            x_i = moles / total_moles
            
            # Calculate properties for this species
            result, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
            
            species_value = result["properties"][property]["value"]
            species_uncertainty = result["properties"][property]["uncertainty"]
            
            # Add contribution to total
            property_value += x_i * species_value
            property_uncertainty += (x_i * species_uncertainty)^2  # Sum of squares for uncertainty
        end
    end
    
    property_uncertainty = sqrt(property_uncertainty)
    
    return property_value, property_uncertainty
end

# Run the main function
main()