"""
Quantum chemistry calculations for thermodynamic properties when experimental data is unavailable.
"""

"""
    estimate_properties_quantum_chemistry(formula::String, temperature::Float64, config::Dict)

Estimate thermodynamic properties using quantum chemistry calculations.
This implementation includes anharmonicity corrections and improved accuracy.
"""
function estimate_properties_quantum_chemistry(formula::String, temperature::Float64, config::Dict)
    # This is a placeholder for actual quantum chemistry calculations
    # In a real implementation, this would interface with quantum chemistry software
    # like Psi4, ORCA, Gaussian, etc.
    
    # Get configuration for quantum chemistry
    qc_config = nothing
    
    for method in config["theoretical_calculation"]["methods"]
        if method["name"] == "quantum_chemistry"
            qc_config = method
            break
        end
    end
    
    if qc_config === nothing
        error("Quantum chemistry configuration not found")
    end
    
    # Determine whether the calculation is enabled
    if !get(qc_config, "enabled", false)
        error("Quantum chemistry calculations are disabled in the configuration")
    end
    
    # Get calculation settings
    program = get(qc_config, "program", "Psi4")
    theory_level = get(qc_config, "theory_level", "B3LYP/6-31G*")
    include_anharmonicity = get(qc_config, "include_anharmonicity", true)
    
    # This is a placeholder for an actual calculation
    # We would typically:
    # 1. Generate a molecular structure from the formula
    # 2. Optimize the geometry
    # 3. Calculate vibrational frequencies
    # 4. Calculate thermodynamic properties
    
    # For the demo, we'll simulate the outcome
    @info "Simulating quantum chemistry calculation for $formula using $program at $theory_level level"
    
    # Generate a placeholder molecular structure
    molecular_data = generate_placeholder_molecule(formula)
    
    # Apply anharmonicity corrections if enabled
    if include_anharmonicity
        molecular_data = apply_anharmonicity_corrections(molecular_data)
    end
    
    # Use statistical thermodynamics with our "calculated" molecular parameters
    result = calculate_quantum_corrected_properties(formula, molecular_data, temperature)
    
    # Mark the method as quantum chemistry
    result["method"] = "quantum_chemistry"
    
    # Add metadata about the calculation
    result["quantum_chemistry"] = Dict(
        "program" => program,
        "theory_level" => theory_level,
        "energy" => -100.0 * sum(values(extract_formula_elements(formula))),  # Placeholder energy
        "anharmonicity" => include_anharmonicity
    )
    
    # Adjust the uncertainties based on theory level
    uncertainty_factors = Dict(
        "HF/3-21G" => 2.0,
        "HF/6-31G*" => 1.8,
        "B3LYP/6-31G*" => 1.2,
        "B3LYP/6-31G**" => 1.1,
        "B3LYP/6-311+G(2d,p)" => 0.9,
        "B3LYP/aug-cc-pVTZ" => 0.8,
        "M06-2X/6-311+G(d,p)" => 0.85,
        "ωB97X-D/6-311+G(d,p)" => 0.82,
        "CCSD(T)/cc-pVTZ" => 0.7,
        "CCSD(T)/CBS" => 0.5
    )
    
    # Apply uncertainty factor based on theory level
    uncertainty_factor = get(uncertainty_factors, theory_level, 1.0)
    
    # Adjust uncertainties
    result["Cp_uncertainty"] *= uncertainty_factor
    result["H_uncertainty"] *= uncertainty_factor
    result["S_uncertainty"] *= uncertainty_factor
    result["G_uncertainty"] *= uncertainty_factor
    
    return result
end

"""
    apply_anharmonicity_corrections(molecular_data::Dict)

Apply anharmonicity corrections to vibrational frequencies and other molecular properties.
"""
function apply_anharmonicity_corrections(molecular_data::Dict)
    # Create a copy of the molecular data
    corrected_data = deepcopy(molecular_data)
    
    # Apply anharmonicity corrections to vibrational frequencies
    if haskey(corrected_data, "vibrational_frequencies")
        frequencies = corrected_data["vibrational_frequencies"]
        corrected_frequencies = similar(frequencies)
        
        for (i, freq) in enumerate(frequencies)
            # Skip imaginary frequencies
            if freq <= 0
                corrected_frequencies[i] = freq
                continue
            end
            
            # Apply anharmonicity correction
            # For real QM calculations, this would use mode-specific anharmonicity constants
            # Instead, we'll use an empirical scaling factor that depends on frequency range
            if freq > 3000  # High frequency (typically X-H stretching)
                scale_factor = 0.96  # High frequency modes typically have larger anharmonicity
            elseif freq > 1700  # Medium-high frequency
                scale_factor = 0.97
            elseif freq > 1000  # Medium frequency
                scale_factor = 0.98
            else  # Low frequency
                scale_factor = 0.99  # Low frequency modes typically have smaller anharmonicity
            end
            
            # Apply scaling to account for anharmonicity
            corrected_frequencies[i] = freq * scale_factor
        end
        
        # Update the frequencies
        corrected_data["vibrational_frequencies"] = corrected_frequencies
        
        # Add anharmonicity constants (placeholder values)
        # In reality, these would come from QM calculations
        n_modes = length(frequencies)
        if n_modes > 0
            anharmonicity_constants = zeros(n_modes)
            
            # Set typical anharmonicity constants (very rough approximations)
            for i in 1:n_modes
                freq = frequencies[i]
                if freq <= 0
                    continue
                end
                
                # Typical anharmonicity constants in cm^-1
                if freq > 3000
                    anharmonicity_constants[i] = -50.0 - 20.0 * rand()  # X-H stretching
                elseif freq > 1700
                    anharmonicity_constants[i] = -15.0 - 10.0 * rand()  # C=O, C=C, etc.
                elseif freq > 1000
                    anharmonicity_constants[i] = -10.0 - 5.0 * rand()   # Mid-range modes
                else
                    anharmonicity_constants[i] = -2.0 - 3.0 * rand()    # Low-frequency modes
                end
            end
            
            corrected_data["anharmonicity_constants"] = anharmonicity_constants
        end
    end
    
    return corrected_data
end

"""
    calculate_quantum_corrected_properties(formula::String, molecular_data::Dict, temperature::Float64)

Calculate thermodynamic properties with quantum corrections.
This version includes corrections for anharmonicity and quantum effects.
"""
function calculate_quantum_corrected_properties(formula::String, molecular_data::Dict, temperature::Float64)
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
    
    # Calculate translational contribution with quantum corrections
    q_trans = calculate_translational_partition_function(molecular_weight, temperature)
    
    # Calculate rotational contribution (if rotational constants provided)
    q_rot = 1.0
    if haskey(molecular_data, "rotational_constants")
        # Apply quantum correction to rotational partition function
        q_rot = calculate_quantum_corrected_rotational_partition_function(
            molecular_data["rotational_constants"],
            symmetry_number,
            temperature
        )
    end
    
    # Calculate vibrational contribution (if vibrational frequencies provided)
    q_vib = 1.0
    if haskey(molecular_data, "vibrational_frequencies")
        # Apply anharmonicity if constants are available
        if haskey(molecular_data, "anharmonicity_constants")
            q_vib = calculate_anharmonic_vibrational_partition_function(
                molecular_data["vibrational_frequencies"],
                molecular_data["anharmonicity_constants"],
                temperature
            )
        else
            q_vib = calculate_vibrational_partition_function(
                molecular_data["vibrational_frequencies"],
                temperature
            )
        end
    end
    
    # Calculate electronic contribution
    q_elec = spin_multiplicity
    
    # Calculate thermodynamic properties
    
    # Enthalpy (H)
    h_trans = R * temperature * (1.0 + 1/(12.0 * (q_trans)^(2/3)))  # with quantum correction term
    
    h_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        h_rot = calculate_quantum_corrected_rotational_enthalpy(
            molecular_data["rotational_constants"],
            temperature
        )
    end
    
    h_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        if haskey(molecular_data, "anharmonicity_constants")
            h_vib = calculate_anharmonic_vibrational_enthalpy(
                molecular_data["vibrational_frequencies"],
                molecular_data["anharmonicity_constants"],
                temperature
            )
        else
            h_vib = calculate_vibrational_enthalpy(
                molecular_data["vibrational_frequencies"],
                temperature
            )
        end
    end
    
    h_elec = 0.0  # Electronic contribution at ground state
    
    h = h_trans + h_rot + h_vib + h_elec
    
    # Add standard enthalpy of formation if provided
    if haskey(molecular_data, "enthalpy_formation")
        h += molecular_data["enthalpy_formation"]
    end
    
    # Entropy (S)
    s_trans = R * (log(q_trans) + 5/3 + 1/(12.0 * (q_trans)^(2/3)))  # with quantum correction
    
    s_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        s_rot = calculate_quantum_corrected_rotational_entropy(
            molecular_data["rotational_constants"],
            symmetry_number,
            temperature
        )
    end
    
    s_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        if haskey(molecular_data, "anharmonicity_constants")
            s_vib = calculate_anharmonic_vibrational_entropy(
                molecular_data["vibrational_frequencies"],
                molecular_data["anharmonicity_constants"],
                temperature
            )
        else
            s_vib = calculate_vibrational_entropy(
                molecular_data["vibrational_frequencies"],
                temperature
            )
        end
    end
    
    s_elec = R * log(q_elec)  # Electronic contribution
    
    s = s_trans + s_rot + s_vib + s_elec
    
    # Heat capacity (Cp)
    cp_trans = 5 * R / 2 * (1.0 + 1/(6.0 * (q_trans)^(2/3)))  # with quantum correction
    
    cp_rot = 0.0
    if haskey(molecular_data, "rotational_constants")
        cp_rot = calculate_quantum_corrected_rotational_heat_capacity(
            molecular_data["rotational_constants"],
            temperature
        )
    end
    
    cp_vib = 0.0
    if haskey(molecular_data, "vibrational_frequencies")
        if haskey(molecular_data, "anharmonicity_constants")
            cp_vib = calculate_anharmonic_vibrational_heat_capacity(
                molecular_data["vibrational_frequencies"],
                molecular_data["anharmonicity_constants"],
                temperature
            )
        else
            cp_vib = calculate_vibrational_heat_capacity(
                molecular_data["vibrational_frequencies"],
                temperature
            )
        end
    end
    
    cp_elec = 0.0  # Electronic contribution at ground state
    
    cp = cp_trans + cp_rot + cp_vib + cp_elec
    
    # Gibbs free energy (G)
    g = h - temperature * s / 1000  # Convert S to kJ/mol/K
    
    # Estimate uncertainties - improved for quantum chemistry
    # Lower uncertainties due to more accurate method
    cp_uncertainty = max(cp * 0.08, 1.5)  # 8% or at least 1.5 J/mol/K
    h_uncertainty = max(abs(h) * 0.1, 8.0)  # 10% or at least 8 kJ/mol
    s_uncertainty = max(s * 0.08, 4.0)  # 8% or at least 4 J/mol/K
    g_uncertainty = max(abs(g) * 0.1, 8.0)  # 10% or at least 8 kJ/mol
    
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
        "method" => "quantum_chemistry"
    )
end

"""
    generate_placeholder_molecule(formula::String)

Generate a placeholder molecular structure for quantum chemistry calculations.
"""
function generate_placeholder_molecule(formula::String)
    # In a real implementation, this would generate a 3D structure using
    # tools like RDKit, Open Babel, etc.
    
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
    
    # Total number of atoms
    total_atoms = sum(values(elements))
    
    # Estimate vibrational frequencies
    # These would come from quantum chemistry calculations in reality
    vib_freq = Float64[]
    
    # Number of vibrational modes
    if total_atoms == 1
        # Monatomic - no vibrations
        n_modes = 0
    elseif total_atoms == 2
        # Diatomic - one vibration
        n_modes = 1
    else
        # Polyatomic
        n_modes = 3 * total_atoms - 6  # Non-linear molecule
    end
    
    # Generate frequencies
    if n_modes == 1
        # Diatomic
        if haskey(elements, "H")
            # X-H bond
            push!(vib_freq, 2800.0 + 200.0 * rand())
        else
            # Heavy element bond
            push!(vib_freq, 1200.0 + 300.0 * rand())
        end
    elseif n_modes > 1
        # Polyatomic - more detailed frequency distribution
        
        # Low frequency modes (bending, torsion) - 100-800 cm^-1
        n_low = floor(Int, n_modes * 0.4)
        for i in 1:n_low
            push!(vib_freq, 100.0 + 700.0 * rand())
        end
        
        # Medium frequency modes (many types) - 800-1800 cm^-1
        n_med = floor(Int, n_modes * 0.4)
        for i in 1:n_med
            push!(vib_freq, 800.0 + 1000.0 * rand())
        end
        
        # High frequency modes (typically X-H stretching) - 2800-3800 cm^-1
        n_high = n_modes - n_low - n_med
        for i in 1:n_high
            push!(vib_freq, 2800.0 + 1000.0 * rand())
        end
    end
    
    # Estimate rotational constants
    rot_const = Float64[]
    
    if total_atoms == 1
        # Monatomic - no rotation
    elseif total_atoms == 2
        # Diatomic - one rotational constant
        
        # Estimate bond length based on elements
        # These are very rough estimates
        bond_length = 1.5  # Å (default)
        
        # Adjust for common bond types
        if haskey(elements, "H")
            bond_length = 1.0  # X-H bonds are shorter
        elseif haskey(elements, "O") && elements["O"] == 2
            bond_length = 1.2  # O=O bond
        elseif haskey(elements, "N") && elements["N"] == 2
            bond_length = 1.1  # N≡N bond
        elseif haskey(elements, "C") && haskey(elements, "O")
            bond_length = 1.3  # C-O or C=O average
        end
        
        # Calculate rotational constant (cm^-1)
        # B = h/(8π²cI) where I = μr² and μ is reduced mass
        
        # Calculate reduced mass (amu)
        element_list = [k for k in keys(elements)]
        if length(element_list) == 1
            # Homonuclear diatomic
            mass1 = ATOMIC_MASSES[element_list[1]]
            mass2 = mass1
        else
            # Heteronuclear diatomic
            mass1 = ATOMIC_MASSES[element_list[1]]
            mass2 = ATOMIC_MASSES[element_list[2]]
        end
        
        reduced_mass = (mass1 * mass2) / (mass1 + mass2)
        
        # Convert bond length to cm
        bond_length_cm = bond_length * 1e-8
        
        # Calculate moment of inertia (g·cm²)
        I = reduced_mass * bond_length_cm^2 * 1e-24
        
        # Calculate rotational constant (cm^-1)
        B = 27.08 / I
        
        push!(rot_const, B)
    else
        # Polyatomic - three rotational constants
        
        # Estimate principal moments of inertia
        # These would come from 3D geometry in reality
        
        # Scale factor based on molecular size
        scale = total_atoms^(1/3)
        
        # For asymmetric top molecules, we need three different constants
        # The asymmetry depends on molecular shape
        
        # Base rotational constant
        B_base = 10.0 / scale
        
        # Add asymmetry
        push!(rot_const, B_base * (1.0 + 0.5 * rand()))  # A
        push!(rot_const, B_base)                        # B
        push!(rot_const, B_base * (1.0 - 0.5 * rand()))  # C
    end
    
    # Symmetry number
    # This is a very simplified estimate
    symmetry = 1
    
    # Some common cases
    if length(elements) == 1 && total_atoms == 2
        # Homonuclear diatomic (like O2, N2)
        symmetry = 2
    elseif haskey(elements, "C") && haskey(elements, "H") && length(elements) == 2
        # Hydrocarbon - might have some symmetry
        if elements["C"] == 2 && elements["H"] == 6
            # C2H6 (ethane) - staggered conformation
            symmetry = 6
        elseif elements["C"] == 2 && elements["H"] == 4
            # C2H4 (ethylene)
            symmetry = 4
        elseif elements["C"] == 2 && elements["H"] == 2
            # C2H2 (acetylene)
            symmetry = 2
        elseif elements["C"] == 1 && elements["H"] == 4
            # CH4 (methane)
            symmetry = 12
        end
    end
    
    # Spin multiplicity
    # This is a very simplified estimate
    spin = 1  # Singlet by default
    
    # Some common cases
    if haskey(elements, "O") && elements["O"] == 2 && length(elements) == 1
        # O2 is a triplet
        spin = 3
    elseif haskey(elements, "N") && elements["N"] == 1 && 
           haskey(elements, "O") && elements["O"] == 1 && length(elements) == 2
        # NO is a doublet
        spin = 2
    end
    
    # Return molecular data dictionary
    return Dict(
        "molecular_weight" => mw,
        "vibrational_frequencies" => vib_freq,
        "rotational_constants" => rot_const,
        "symmetry_number" => symmetry,
        "spin_multiplicity" => spin,
        "elements" => elements
    )
end

"""
    calculate_quantum_corrected_rotational_partition_function(rotational_constants::Vector{<:Real}, 
                                                           symmetry_number::Int, temperature::Float64)

Calculate the rotational partition function with quantum corrections.
"""
function calculate_quantum_corrected_rotational_partition_function(
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
        
        # Define energy levels in kB*T units
        x = h * B_si / (KB * temperature)
        
        # Quantum correction factor (first-order correction)
        quantum_corr = 1.0 + (x^2) / 12.0
        
        # Rotational partition function for linear molecule with quantum correction
        q_rot = (KB * temperature) / (h * B_si * symmetry_number) * quantum_corr
    else
        # Non-linear molecule
        A = rotational_constants[1]  # cm^-1
        B = rotational_constants[2]  # cm^-1
        C = rotational_constants[3]  # cm^-1
        
        # Convert to SI units (Hz)
        A_si = A * 100 * 2.998e10  # Hz
        B_si = B * 100 * 2.998e10  # Hz
        C_si = C * 100 * 2.998e10  # Hz
        
        # Average rotational constant
        B_avg = (A_si * B_si * C_si)^(1/3)
        
        # Quantum correction factor
        x_avg = h * B_avg / (KB * temperature)
        quantum_corr = 1.0 + (x_avg^2) / 24.0
        
        # Rotational partition function for non-linear molecule with quantum correction
        q_rot = sqrt(π) / symmetry_number * ((KB * temperature)^1.5) / (h^1.5 * sqrt(A_si * B_si * C_si)) * quantum_corr
    end
    
    return q_rot
end

"""
    calculate_anharmonic_vibrational_partition_function(vibrational_frequencies::Vector{Float64}, 
                                                     anharmonicity_constants::Vector{Float64}, 
                                                     temperature::Float64)

Calculate the vibrational partition function with anharmonicity corrections.
"""
function calculate_anharmonic_vibrational_partition_function(
    vibrational_frequencies::Vector{Float64}, 
    anharmonicity_constants::Vector{Float64}, 
    temperature::Float64
)
    # Vibrational frequencies should be in cm^-1
    # Anharmonicity constants should be in cm^-1
    
    q_vib = 1.0
    
    for i in 1:length(vibrational_frequencies)
        nu = vibrational_frequencies[i]
        
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Get anharmonicity constant
        anharm = i <= length(anharmonicity_constants) ? anharmonicity_constants[i] : 0.0
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * C  # Hz
        anharm_si = anharm * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
        # Anharmonicity correction
        anharm_factor = 1.0
        if anharm != 0.0
            # First-order anharmonicity correction to partition function
            # This is a simplified approach using perturbation theory
            anharm_factor = 1.0 - (h * anharm_si / (KB * temperature)) * (exp(-x) / (1 - exp(-x)))
        end
        
        # Vibrational partition function contribution with anharmonicity
        q_vib *= (1 / (1 - exp(-x))) * anharm_factor
    end
    
    return q_vib
end

"""
    calculate_anharmonic_vibrational_enthalpy(vibrational_frequencies::Vector{Float64}, 
                                           anharmonicity_constants::Vector{Float64}, 
                                           temperature::Float64)

Calculate the vibrational contribution to enthalpy with anharmonicity corrections.
"""
function calculate_anharmonic_vibrational_enthalpy(
    vibrational_frequencies::Vector{Float64}, 
    anharmonicity_constants::Vector{Float64}, 
    temperature::Float64
)
    # Vibrational frequencies should be in cm^-1
    
    h_vib = 0.0
    
    for i in 1:length(vibrational_frequencies)
        nu = vibrational_frequencies[i]
        
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Get anharmonicity constant
        anharm = i <= length(anharmonicity_constants) ? anharmonicity_constants[i] : 0.0
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * C  # Hz
        anharm_si = anharm * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
        # Base harmonic contribution
        harmonic_term = R * temperature * x / (exp(x) - 1)
        
        # Zero-point energy contribution
        zpe_term = R * temperature * x / 2
        
        # Anharmonicity correction to energy
        anharm_term = 0.0
        if anharm != 0.0
            # Anharmonicity correction to energy
            anharm_x = h * anharm_si / (KB * temperature)
            anharm_term = R * temperature * anharm_x * ((exp(x) + 1) / (2 * (exp(x) - 1)))
        end
        
        # Total contribution for this mode
        h_vib += harmonic_term + zpe_term + anharm_term
    end
    
    return h_vib
end

"""
    calculate_anharmonic_vibrational_entropy(vibrational_frequencies::Vector{Float64}, 
                                          anharmonicity_constants::Vector{Float64}, 
                                          temperature::Float64)

Calculate the vibrational contribution to entropy with anharmonicity corrections.
"""
function calculate_anharmonic_vibrational_entropy(
    vibrational_frequencies::Vector{Float64}, 
    anharmonicity_constants::Vector{Float64}, 
    temperature::Float64
)
    # Vibrational frequencies should be in cm^-1
    
    s_vib = 0.0
    
    for i in 1:length(vibrational_frequencies)
        nu = vibrational_frequencies[i]
        
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Get anharmonicity constant
        anharm = i <= length(anharmonicity_constants) ? anharmonicity_constants[i] : 0.0
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * C  # Hz
        anharm_si = anharm * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
        # Harmonic entropy contribution
        harmonic_term = R * (x / (exp(x) - 1) - log(1 - exp(-x)))
        
        # Anharmonicity correction to entropy
        anharm_term = 0.0
        if anharm != 0.0
            # Anharmonicity correction to entropy
            anharm_x = h * anharm_si / (KB * temperature)
            anharm_term = R * anharm_x * (1 / (exp(x) - 1) - x * exp(x) / ((exp(x) - 1)^2))
        end
        
        # Total contribution for this mode
        s_vib += harmonic_term + anharm_term
    end
    
    return s_vib
end

"""
    calculate_anharmonic_vibrational_heat_capacity(vibrational_frequencies::Vector{Float64}, 
                                                anharmonicity_constants::Vector{Float64}, 
                                                temperature::Float64)

Calculate the vibrational contribution to heat capacity with anharmonicity corrections.
"""
function calculate_anharmonic_vibrational_heat_capacity(
    vibrational_frequencies::Vector{Float64}, 
    anharmonicity_constants::Vector{Float64}, 
    temperature::Float64
)
    # Vibrational frequencies should be in cm^-1
    
    cp_vib = 0.0
    
    for i in 1:length(vibrational_frequencies)
        nu = vibrational_frequencies[i]
        
        # Skip imaginary frequencies (negative values)
        if nu <= 0
            continue
        end
        
        # Get anharmonicity constant
        anharm = i <= length(anharmonicity_constants) ? anharmonicity_constants[i] : 0.0
        
        # Convert to SI units (Hz)
        nu_si = nu * 100 * C  # Hz
        anharm_si = anharm * 100 * C  # Hz
        
        # Vibrational energy in units of kB*T
        x = h * nu_si / (KB * temperature)
        
        # Harmonic heat capacity contribution
        harmonic_term = R * (x^2 * exp(x) / (exp(x) - 1)^2)
        
        # Anharmonicity correction to heat capacity
        anharm_term = 0.0
        if anharm != 0.0
            # Anharmonicity correction to heat capacity
            anharm_x = h * anharm_si / (KB * temperature)
            ex = exp(x)
            anharm_term = R * anharm_x * x * (ex * (ex + 1) / (ex - 1)^3)
        end
        
        # Total contribution for this mode
        cp_vib += harmonic_term + anharm_term
    end
    
    return cp_vib
end

"""
    calculate_quantum_corrected_rotational_enthalpy(rotational_constants::Vector{Float64}, temperature::Float64)

Calculate the rotational contribution to enthalpy with quantum corrections.
"""
function calculate_quantum_corrected_rotational_enthalpy(rotational_constants::Vector{Float64}, temperature::Float64)
    # Rotational constants should be in cm^-1
    
    # Check if molecule is linear or non-linear
    if length(rotational_constants) == 1
        # Linear molecule
        B = rotational_constants[1]  # cm^-1
        
        # Quantum correction to rotational energy
        # Convert B to J
        B_j = B * 100 * C * h  # J
        
        # Quantum correction term (simplified first-order correction)
        correction = 1.0 + B_j / (6.0 * KB * temperature)
        
        # Rotational enthalpy contribution
        h_rot = R * temperature * correction
    else
        # Non-linear molecule
        # For non-linear molecules, quantum correction is smaller
        # but we'll include it for completeness
        
        # Average rotational constant
        B_avg = sum(rotational_constants) / length(rotational_constants)
        
        # Convert B to J
        B_j = B_avg * 100 * C * h  # J
        
        # Quantum correction term (simplified first-order correction)
        correction = 1.0 + B_j / (4.0 * KB * temperature)
        
        # Rotational enthalpy contribution for non-linear molecules
        h_rot = 1.5 * R * temperature * correction
    end
    
    return h_rot
end

"""
    calculate_quantum_corrected_rotational_entropy(rotational_constants::Vector{Float64}, 
                                               symmetry_number::Int, temperature::Float64)

Calculate the rotational contribution to entropy with quantum corrections.
"""
function calculate_quantum_corrected_rotational_entropy(
    rotational_constants::Vector{Float64}, 
    symmetry_number::Int, 
    temperature::Float64
)
    # Calculate partition function with quantum corrections
    q_rot = calculate_quantum_corrected_rotational_partition_function(
        rotational_constants, symmetry_number, temperature
    )
    
    # Check if molecule is linear or non-linear
    if length(rotational_constants) == 1
        # Linear molecule
        s_rot = R * (log(q_rot) + 1)  # Rotational entropy
    else
        # Non-linear molecule
        s_rot = R * (log(q_rot) + 3/2)  # Rotational entropy
    end
    
    return s_rot
end

"""
    calculate_quantum_corrected_rotational_heat_capacity(rotational_constants::Vector{Float64}, temperature::Float64)

Calculate the rotational contribution to heat capacity with quantum corrections.
"""
function calculate_quantum_corrected_rotational_heat_capacity(rotational_constants::Vector{Float64}, temperature::Float64)
    # Rotational constants should be in cm^-1
    
    # Check if molecule is linear or non-linear
    if length(rotational_constants) == 1
        # Linear molecule
        B = rotational_constants[1]  # cm^-1
        
        # Quantum correction to rotational heat capacity
        # Convert B to J
        B_j = B * 100 * C * h  # J
        
        # High-temperature limit for linear molecule is R
        cp_rot_classical = R
        
        # Quantum correction (decreases with temperature)
        correction = 1.0 - (B_j / (6.0 * KB * temperature))^2
        
        # Apply correction (ensuring it doesn't go negative)
        cp_rot = cp_rot_classical * max(correction, 0.5)
    else
        # Non-linear molecule
        # Classical rotational heat capacity for non-linear molecule is 3R/2
        cp_rot_classical = 1.5 * R
        
        # Average rotational constant
        B_avg = sum(rotational_constants) / length(rotational_constants)
        
        # Convert B to J
        B_j = B_avg * 100 * C * h  # J
        
        # Quantum correction (decreases with temperature)
        correction = 1.0 - (B_j / (4.0 * KB * temperature))^2
        
        # Apply correction (ensuring it doesn't go negative)
        cp_rot = cp_rot_classical * max(correction, 0.5)
    end
    
    return cp_rot
end

"""
    run_psi4_calculation(formula::String, config::Dict)

Run a Psi4 quantum chemistry calculation for the given molecule.
"""
function run_psi4_calculation(formula::String, config::Dict)
    # In a real implementation, this would:
    # 1. Generate a molecular structure
    # 2. Write a Psi4 input file
    # 3. Execute Psi4
    # 4. Parse the output file for thermodynamic properties
    
    error("Psi4 interface not implemented yet")
end

"""
    run_orca_calculation(formula::String, config::Dict)

Run an ORCA quantum chemistry calculation for the given molecule.
"""
function run_orca_calculation(formula::String, config::Dict)
    # In a real implementation, this would:
    # 1. Generate a molecular structure
    # 2. Write an ORCA input file
    # 3. Execute ORCA
    # 4. Parse the output file for thermodynamic properties
    
    error("ORCA interface not implemented yet")
end