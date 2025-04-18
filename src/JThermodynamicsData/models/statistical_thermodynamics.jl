"""
Statistical thermodynamics methods for calculating thermodynamic properties from molecular data.
"""

# Add lowercase constant aliases to match usage in functions
const h = H
const kb = KB
const na = NA
const c = C

"""
    estimate_properties_statistical_thermodynamics(formula::String, molecular_data::Dict, temperature::Float64)

Estimate thermodynamic properties using statistical thermodynamics.
"""
function estimate_properties_statistical_thermodynamics(
    formula::String, 
    molecular_data::Dict, 
    temperature::Float64
)
    # Statistical thermodynamics calculates thermodynamic properties from
    # molecular properties like vibrational frequencies, rotational constants, etc.
    
    # Check if we have the required molecular data
    required_keys = ["molecular_weight", "symmetry_number", "spin_multiplicity"]
    
    for key in required_keys
        if !haskey(molecular_data, key)
            error("Missing required molecular data: $key")
        end
    end
    
    # Get molecular properties
    molecular_weight = molecular_data["molecular_weight"]  # g/mol
    symmetry_number = molecular_data["symmetry_number"]
    spin_multiplicity = molecular_data["spin_multiplicity"]
    
    # Calculate translational contribution
    q_trans = calculate_translational_partition_function(molecular_weight, temperature)
    
    # Calculate rotational contribution (if rotational constants provided)
    q_rot = 1.0
    if haskey(molecular_data, "rotational_constants")
        q_rot = calculate_rotational_partition_function(
            molecular_data["rotational_constants"],
            symmetry_number,
            temperature
        )
    end
    
    # Calculate vibrational contribution (if vibrational frequencies provided)
    q_vib = 1.0
    if haskey(molecular_data, "vibrational_frequencies")
        q_vib = calculate_vibrational_partition_function(
            molecular_data["vibrational_frequencies"],
            temperature
        )
    end
    
    # Calculate electronic contribution
    q_elec = spin_multiplicity
    
    # Calculate thermodynamic properties
    
    # Enthalpy (H)
    h_trans = R * temperature  # Translational contribution to enthalpy
    
    h_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        h_rot = R * temperature  # Rotational contribution for non-linear molecules
    end
    
    h_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        h_vib = calculate_vibrational_enthalpy(
            molecular_data["vibrational_frequencies"],
            temperature
        )
    end
    
    h_elec = 0.0  # Electronic contribution at ground state
    
    h = h_trans + h_rot + h_vib + h_elec
    
    # Add standard enthalpy of formation if provided
    if haskey(molecular_data, "enthalpy_formation")
        h += molecular_data["enthalpy_formation"]
    end
    
    # Entropy (S)
    s_trans = R * (log(q_trans) + 1 + 1.5)  # Translational contribution to entropy
    
    s_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        s_rot = R * (log(q_rot) + 1)  # Rotational contribution
    end
    
    s_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        s_vib = calculate_vibrational_entropy(
            molecular_data["vibrational_frequencies"],
            temperature
        )
    end
    
    s_elec = R * log(q_elec)  # Electronic contribution
    
    s = s_trans + s_rot + s_vib + s_elec
    
    # Heat capacity (Cp)
    cp_trans = 5 * R / 2  # Translational contribution to heat capacity
    
    cp_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        cp_rot = R  # Rotational contribution for non-linear molecules (3/2 * R for linear)
    end
    
    cp_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        cp_vib = calculate_vibrational_heat_capacity(
            molecular_data["vibrational_frequencies"],
            temperature
        )
    end
    
    cp_elec = 0.0  # Electronic contribution at ground state
    
    cp = cp_trans + cp_rot + cp_vib + cp_elec
    
    # Gibbs free energy (G)
    g = h - temperature * s / 1000  # Convert S to kJ/mol/K
    
    # Estimate uncertainties
    # Statistical thermodynamics uncertainties depend on the accuracy of the molecular data
    cp_uncertainty = max(cp * 0.1, 2.0)  # 10% or at least 2 J/mol/K
    h_uncertainty = max(abs(h) * 0.15, 10.0)  # 15% or at least 10 kJ/mol
    s_uncertainty = max(s * 0.12, 5.0)  # 12% or at least 5 J/mol/K
    g_uncertainty = max(abs(g) * 0.15, 10.0)  # 15% or at least 10 kJ/mol
    
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
        "method" => "statistical_thermodynamics"
    )
end

"""
    calculate_translational_partition_function(molecular_weight::Float64, temperature::Float64)

Calculate the translational partition function.
"""
function calculate_translational_partition_function(molecular_weight::Float64, temperature::Float64)
    # Convert molecular weight from g/mol to kg/molecule
    m = molecular_weight / (NA * 1000)
    
    # Thermal wavelength
    lambda = h / sqrt(2 * π * m * KB * temperature)
    
    # Volume per molecule (approximation for ideal gas at standard conditions)
    V = (R * temperature) / (P_STANDARD * NA)
    
    # Translational partition function
    q_trans = V / lambda^3
    
    return q_trans
end

"""
    calculate_rotational_partition_function(rotational_constants::Vector{Float64}, 
                                          symmetry_number::Int, temperature::Float64)

Calculate the rotational partition function.
"""
function calculate_rotational_partition_function(
    rotational_constants::Vector{<:Real}, 
    symmetry_number::Int, 
    temperature::Float64
)
    # Rotational constants should be in cm^-1
    
    # Check if molecule is linear or non-linear
    if length(rotational_constants) == 1
        # Linear molecule
        B = rotational_constants[1]  # cm^-1
        
        # Convert to SI units (Hz)
        B_si = B * C * 100  # Hz
        
        # Rotational partition function for linear molecule
        q_rot = (KB * temperature) / (h * B_si * symmetry_number)
    else
        # Non-linear molecule
        A = rotational_constants[1]  # cm^-1
        B = rotational_constants[2]  # cm^-1
        C = rotational_constants[3]  # cm^-1
        
        # Convert to SI units (Hz)
        A_si = A * 100 * 2.998e10  # Hz
        B_si = B * 100 * 2.998e10  # Hz
        C_si = C * 100 * 2.998e10  # Hz
        
        # Rotational partition function for non-linear molecule
        q_rot = sqrt(π) / symmetry_number * ((KB * temperature)^1.5) / (h^1.5 * sqrt(A_si * B_si * C_si))
    end
    
    return q_rot
end

"""
    calculate_vibrational_partition_function(vibrational_frequencies::Vector{Float64}, temperature::Float64)

Calculate the vibrational partition function.
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
        nu_si = nu * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
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
        nu_si = nu * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
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
        nu_si = nu * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
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
        nu_si = nu * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
        # Vibrational heat capacity contribution
        cp_vib += R * (x^2 * exp(x) / (exp(x) - 1)^2)
    end
    
    return cp_vib
end

"""
    estimate_molecular_properties(formula::String)

Estimate basic molecular properties from the chemical formula.
For use when experimental data isn't available.
"""
function estimate_molecular_properties(formula::String)
    # Parse formula
    elements = extract_formula_elements(formula)
    
    # Estimate molecular weight
    mw = 0.0
    for (element, count) in elements
        if haskey(ATOMIC_MASSES, element)
            mw += ATOMIC_MASSES[element] * count
        else
            @warn "Unknown element in formula: $element"
        end
    end
    
    # Very rough estimation of molecular properties
    # In practice, these would come from quantum chemistry calculations or experiments
    
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
        append!(vib_freq, rand(400.0:800.0, n_low))
        
        # Medium frequency modes (800-1500 cm^-1)
        n_medium = floor(Int, n_modes / 3)
        append!(vib_freq, rand(800.0:1500.0, n_medium))
        
        # High frequency modes (1500-3500 cm^-1)
        n_high = n_modes - n_low - n_medium
        append!(vib_freq, rand(1500.0:3500.0, n_high))
        
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
        B = h / (8 * π^2 * C * I) / 100  # cm^-1
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
        B = h / (8 * π^2 * C * I) / 100  # cm^-1
        
        # For non-spherical molecules, introduce some asymmetry
        push!(rot_const, B * (1.0 + 0.2 * rand()))
        push!(rot_const, B * (1.0 - 0.1 * rand()))
        push!(rot_const, B * (1.0 - 0.3 * rand()))
    end
    
    # Return estimated molecular properties
    return Dict(
        "molecular_weight" => mw,
        "symmetry_number" => symmetry,
        "spin_multiplicity" => spin,
        "vibrational_frequencies" => vib_freq,
        "rotational_constants" => rot_const
    )
end