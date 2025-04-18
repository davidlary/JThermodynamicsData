#!/usr/bin/env julia

"""
Reaction Simulator for Thermodynamic Calculations

This script provides tools for simulating chemical reactions based on the hierarchical thermodynamics calculator.
It includes functionality for reaction equilibrium, reaction kinetics, and reactive mixture properties.
"""

module ReactionSimulator

using Plots
using Printf
using Statistics
using LinearAlgebra
using Optim

# Import hierarchical calculator for thermodynamic data
include("hierarchical_calculator.jl")

# Universal gas constant
const R = 8.31446261815324  # J/(mol·K)

export parse_reaction_equation, calculate_delta_g_reaction, calculate_equilibrium_constant,
       calculate_reaction_equilibrium, calculate_mixed_reactant_properties, calculate_arrhenius_rate,
       simulate_batch_reaction, simulate_temperature_dependence, plot_reaction_equilibrium

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
    
    # Function to parse one side of the equation
    function parse_side(side_str::String)
        species = Dict{String, Float64}()
        
        # Split by + sign and trim whitespace
        components = split(side_str, "+")
        
        for comp in components
            comp = strip(comp)
            if isempty(comp)
                continue
            end
            
            # Extract stoichiometric coefficient and species name
            m = match(r"^(\d*\.?\d*)\s*([A-Za-z0-9\(\)]+)$", comp)
            
            if m === nothing
                m = match(r"^([A-Za-z0-9\(\)]+)$", comp)
                if m === nothing
                    error("Invalid component in reaction: $comp")
                else
                    # No coefficient specified, default to 1
                    species_name = m[1]
                    coef = 1.0
                end
            else
                coef_str = m[1]
                species_name = m[2]
                
                # Default coefficient is 1 if not specified
                coef = isempty(coef_str) ? 1.0 : parse(Float64, coef_str)
            end
            
            # Add to the dictionary
            if haskey(species, species_name)
                species[species_name] += coef
            else
                species[species_name] = coef
            end
        end
        
        return species
    end
    
    reactants = parse_side(reactants_str)
    products = parse_side(products_str)
    
    return reactants, products
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
        result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
        g_value = result["properties"]["G"]["value"]
        g_uncertainty = result["properties"]["G"]["uncertainty"]
        
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
        result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
        g_value = result["properties"]["G"]["value"]
        g_uncertainty = result["properties"]["G"]["uncertainty"]
        
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
    calculate_delta_h_reaction(reactants::Dict{String, Float64}, products::Dict{String, Float64}, 
                             temperature::Float64, data_sources::Dict)

Calculate the enthalpy change for a reaction.
"""
function calculate_delta_h_reaction(reactants::Dict{String, Float64}, products::Dict{String, Float64}, 
                                  temperature::Float64, data_sources::Dict)
    # Calculate enthalpy for reactants
    h_reactants = 0.0
    h_reactants_uncertainty = 0.0
    
    for (species, coef) in reactants
        # Get formula (assuming species name is formula for simplicity)
        formula = species
        
        # Calculate properties for this species
        result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
        h_value = result["properties"]["H"]["value"]
        h_uncertainty = result["properties"]["H"]["uncertainty"]
        
        # Add contribution to total
        h_reactants += coef * h_value
        h_reactants_uncertainty += (coef * h_uncertainty)^2  # Sum of squares for uncertainty
    end
    h_reactants_uncertainty = sqrt(h_reactants_uncertainty)
    
    # Calculate enthalpy for products
    h_products = 0.0
    h_products_uncertainty = 0.0
    
    for (species, coef) in products
        # Get formula (assuming species name is formula for simplicity)
        formula = species
        
        # Calculate properties for this species
        result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
        h_value = result["properties"]["H"]["value"]
        h_uncertainty = result["properties"]["H"]["uncertainty"]
        
        # Add contribution to total
        h_products += coef * h_value
        h_products_uncertainty += (coef * h_uncertainty)^2  # Sum of squares for uncertainty
    end
    h_products_uncertainty = sqrt(h_products_uncertainty)
    
    # Calculate ΔH = H(products) - H(reactants)
    delta_h = h_products - h_reactants
    delta_h_uncertainty = sqrt(h_products_uncertainty^2 + h_reactants_uncertainty^2)
    
    return delta_h, delta_h_uncertainty
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
                result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
                g_std = result["properties"]["G"]["value"]
                
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
            result, _, _ = progressively_refine_thermodynamic_data(species, formula, temperature, data_sources)
            
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

"""
    calculate_arrhenius_rate(a::Float64, e_act::Float64, temperature::Float64)

Calculate the reaction rate constant using the Arrhenius equation.
"""
function calculate_arrhenius_rate(a::Float64, e_act::Float64, temperature::Float64)
    # Arrhenius equation: k = A * exp(-Ea/(RT))
    # A: pre-exponential factor
    # Ea: activation energy in kJ/mol
    # R: gas constant in kJ/(mol·K)
    # T: temperature in K
    
    # Convert R to kJ/(mol·K) for consistency with e_act in kJ/mol
    r_kj = R / 1000.0  # 0.008314 kJ/(mol·K)
    
    # Calculate rate constant
    k = a * exp(-e_act / (r_kj * temperature))
    
    return k
end

"""
    simulate_batch_reaction(equation::String, initial_composition::Dict{String, Float64}, 
                          time_span::Vector{Float64}, temperature::Float64, 
                          a_forward::Float64, e_act_forward::Float64)

Simulate a batch reaction using kinetic rate equations.
"""
function simulate_batch_reaction(equation::String, initial_composition::Dict{String, Float64}, 
                               time_span::Vector{Float64}, temperature::Float64, 
                               a_forward::Float64, e_act_forward::Float64)
    # Parse the reaction equation
    reactants, products = parse_reaction_equation(equation)
    
    # Calculate the rate constant
    k_forward = calculate_arrhenius_rate(a_forward, e_act_forward, temperature)
    
    # All species in the reaction
    all_species = Set{String}()
    for species in keys(reactants)
        push!(all_species, species)
    end
    for species in keys(products)
        push!(all_species, species)
    end
    
    # Ensure all species have initial moles defined
    composition = Dict{String, Float64}()
    for species in all_species
        if haskey(initial_composition, species)
            composition[species] = initial_composition[species]
        else
            composition[species] = 0.0
        end
    end
    
    # Create arrays to store results
    times = Float64[]
    compositions = Dict{String, Vector{Float64}}()
    for species in all_species
        compositions[species] = Float64[]
    end
    
    # Function to calculate reaction rate
    function reaction_rate(composition::Dict{String, Float64})
        # Calculate rate based on reactant concentrations
        rate = k_forward
        
        for (species, coef) in reactants
            # Rate depends on concentration of each reactant raised to its stoichiometric coefficient
            rate *= composition[species]^coef
        end
        
        return rate
    end
    
    # Time step
    dt = (time_span[2] - time_span[1]) / 1000
    
    # Simulate the reaction
    t = time_span[1]
    push!(times, t)
    for species in all_species
        push!(compositions[species], composition[species])
    end
    
    while t < time_span[2]
        # Calculate current reaction rate
        rate = reaction_rate(composition)
        
        # Update all species
        for species in all_species
            # Change in moles = ±rate * dt * stoichiometric coefficient
            if haskey(reactants, species)
                composition[species] -= rate * dt * reactants[species]
            end
            if haskey(products, species)
                composition[species] += rate * dt * products[species]
            end
            
            # Ensure non-negative concentrations
            composition[species] = max(0.0, composition[species])
        end
        
        # Update time
        t += dt
        push!(times, t)
        for species in all_species
            push!(compositions[species], composition[species])
        end
    end
    
    return times, compositions
end

"""
    simulate_temperature_dependence(equation::String, initial_composition::Dict{String, Float64},
                                  temp_range::Vector{Float64}, temp_step::Float64, pressure::Float64,
                                  data_sources::Dict)

Simulate the effect of temperature on reaction equilibrium.
"""
function simulate_temperature_dependence(equation::String, initial_composition::Dict{String, Float64},
                                       temp_range::Vector{Float64}, temp_step::Float64, pressure::Float64,
                                       data_sources::Dict)
    # Parse the reaction equation
    reactants, products = parse_reaction_equation(equation)
    
    # Generate temperature points
    temperatures = temp_range[1]:temp_step:temp_range[2]
    
    # Initialize arrays for results
    delta_g_values = Float64[]
    delta_h_values = Float64[]
    k_values = Float64[]
    extents = Float64[]
    
    # Compositions at each temperature
    compositions = Dict{Float64, Dict{String, Float64}}()
    
    for temp in temperatures
        # Calculate thermodynamic properties
        delta_g, _ = calculate_delta_g_reaction(reactants, products, temp, data_sources)
        delta_h, _ = calculate_delta_h_reaction(reactants, products, temp, data_sources)
        k = calculate_equilibrium_constant(delta_g, temp)
        
        # Calculate equilibrium
        equilibrium_comp, extent = calculate_reaction_equilibrium(
            equation, temp, pressure, initial_composition, data_sources
        )
        
        # Store results
        push!(delta_g_values, delta_g)
        push!(delta_h_values, delta_h)
        push!(k_values, k)
        push!(extents, extent)
        compositions[temp] = equilibrium_comp
    end
    
    return temperatures, delta_g_values, delta_h_values, k_values, extents, compositions
end

"""
    plot_reaction_equilibrium(temperatures, delta_g_values, delta_h_values, k_values, extents, 
                            equation::String, use_log_scale::Bool=true)

Create plots for reaction equilibrium data.
"""
function plot_reaction_equilibrium(temperatures, delta_g_values, delta_h_values, k_values, extents, 
                                 equation::String, use_log_scale::Bool=true)
    # Plot Gibbs energy change
    p1 = plot(
        temperatures,
        delta_g_values,
        xlabel="Temperature (K)",
        ylabel="ΔG° (kJ/mol)",
        label="ΔG°",
        title="Gibbs Energy Change",
        lw=2,
        grid=true,
        legend=:topright
    )
    
    # Plot enthalpy change
    p2 = plot(
        temperatures,
        delta_h_values,
        xlabel="Temperature (K)",
        ylabel="ΔH° (kJ/mol)",
        label="ΔH°",
        title="Enthalpy Change",
        lw=2,
        grid=true,
        legend=:topright
    )
    
    # Plot equilibrium constant (log scale)
    p3 = plot(
        temperatures,
        log10.(k_values),
        xlabel="Temperature (K)",
        ylabel="log₁₀(K)",
        label="log₁₀(K)",
        title="Equilibrium Constant",
        lw=2,
        grid=true,
        legend=:topright
    )
    
    # Plot reaction extent
    p4 = plot(
        temperatures,
        extents,
        xlabel="Temperature (K)",
        ylabel="Reaction Extent",
        label="Extent",
        title="Reaction Extent for\n$equation",
        lw=2,
        grid=true,
        legend=:topright
    )
    
    # Apply logarithmic scale to temperature axis if requested
    if use_log_scale
        p1 = plot!(p1, xscale=:log10)
        p2 = plot!(p2, xscale=:log10)
        p3 = plot!(p3, xscale=:log10)
        p4 = plot!(p4, xscale=:log10)
    end
    
    # Combine all plots
    p = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800), dpi=300)
    
    return p
end

# Make hierarchical calculator functions available in this module
import Main.progressively_refine_thermodynamic_data
import Main.calculate_properties_range
import Main.fit_nasa7_coefficients
import Main.calculate_theoretical_properties
import Main.calculate_nasa_properties
import Main.load_polynomial_data

end # module ReactionSimulator

# If this script is run directly, provide an example
if abspath(PROGRAM_FILE) == @__FILE__
    using .ReactionSimulator
    using Plots
    using Printf
    
    println("Reaction Simulator for Thermodynamic Calculations")
    println("=================================================")
    
    # Load thermodynamic data
    data_sources = load_polynomial_data()
    
    # Example: Methane combustion: CH4 + 2O2 = CO2 + 2H2O
    reaction_equation = "CH4 + 2O2 = CO2 + 2H2O"
    temperature = 1000.0  # K
    pressure = 1.0  # bar
    
    # Initial composition (moles)
    initial_composition = Dict(
        "CH4" => 1.0,
        "O2" => 2.0,
        "CO2" => 0.0,
        "H2O" => 0.0
    )
    
    println("\nExample: Methane combustion reaction")
    println("Reaction: $reaction_equation")
    println("Temperature: $temperature K")
    println("Pressure: $pressure bar")
    
    # Calculate standard Gibbs energy change
    reactants, products = parse_reaction_equation(reaction_equation)
    delta_g, delta_g_uncertainty = calculate_delta_g_reaction(reactants, products, temperature, data_sources)
    delta_h, delta_h_uncertainty = calculate_delta_h_reaction(reactants, products, temperature, data_sources)
    
    # Calculate equilibrium constant
    K = calculate_equilibrium_constant(delta_g, temperature)
    
    println("\nThermodynamic analysis:")
    @printf("  ΔG° = %.2f ± %.2f kJ/mol\n", delta_g, delta_g_uncertainty)
    @printf("  ΔH° = %.2f ± %.2f kJ/mol\n", delta_h, delta_h_uncertainty)
    @printf("  K = %.4e\n", K)
    
    # Calculate equilibrium composition
    equil_composition, extent = calculate_reaction_equilibrium(
        reaction_equation, temperature, pressure, initial_composition, data_sources
    )
    
    println("\nEquilibrium composition (moles):")
    for (species, moles) in sort(collect(equil_composition), by=x->x[1])
        @printf("  %s: %.6f\n", species, moles)
    end
    @printf("  Reaction extent: %.6f\n", extent)
    
    # Simulate temperature dependence
    println("\nSimulating temperature dependence...")
    temps, delta_g_values, delta_h_values, k_values, extents, compositions = 
        simulate_temperature_dependence(
            reaction_equation, initial_composition, [500.0, 1500.0], 100.0, pressure, data_sources
        )
    
    # Create equilibrium plots
    p = plot_reaction_equilibrium(
        temps, delta_g_values, delta_h_values, k_values, extents, reaction_equation
    )
    
    # Create output directory if it doesn't exist
    plot_dir = joinpath(@__DIR__, "plots")
    mkpath(plot_dir)
    
    # Save the plot
    plot_file = joinpath(plot_dir, "methane_combustion_equilibrium.png")
    savefig(p, plot_file)
    println("  Plot saved to: $plot_file")
    
    # Simulate reaction kinetics
    println("\nSimulating reaction kinetics...")
    a_forward = 1.0e10  # Arbitrary pre-exponential factor
    e_act_forward = 100.0  # Arbitrary activation energy in kJ/mol
    
    time_span = [0.0, 10.0]  # seconds
    times, compositions = simulate_batch_reaction(
        reaction_equation, initial_composition, time_span, temperature, a_forward, e_act_forward
    )
    
    # Create kinetics plot
    p_kinetics = plot(
        title="Methane Combustion Kinetics at $temperature K",
        xlabel="Time (s)",
        ylabel="Moles",
        grid=true,
        size=(800, 600),
        dpi=300
    )
    
    # Add each species to the plot
    for (species, moles) in compositions
        plot!(p_kinetics, times, moles, label=species, lw=2)
    end
    
    # Save the kinetics plot
    kinetics_file = joinpath(plot_dir, "methane_combustion_kinetics.png")
    savefig(p_kinetics, kinetics_file)
    println("  Kinetics plot saved to: $kinetics_file")
    
    println("\nReaction simulation complete!")
end