# Hierarchical Thermodynamics Calculator

This package provides a standalone implementation of the JThermodynamics hierarchical approach to thermodynamic data calculation, combined with a reaction simulator for chemical equilibrium and kinetics calculations.

## Introduction

The hierarchical approach to thermodynamic data starts with theoretical estimates and then progressively refines through increasingly accurate experimental sources. This implementation addresses several key needs:

1. Provides a complete, standalone calculator that doesn't depend on database infrastructure
2. Maintains proper uncertainty propagation throughout calculations
3. Produces NASA-7 polynomials for both low and high temperature ranges
4. Includes visualization tools for comparing different calculation methods
5. Adds reaction equilibrium and kinetics capabilities

## Features

### Hierarchical Thermodynamic Data Calculation

- Theoretical calculations using statistical thermodynamics and group contribution methods
- Progressive refinement through a hierarchy of data sources with different priorities
- Proper uncertainty estimation and propagation
- NASA-7 polynomial fitting from calculated data
- Visualization of thermodynamic properties with uncertainty bands

### Reaction Simulation

- Parsing and validation of chemical reaction equations
- Calculation of reaction thermodynamics (ΔG, ΔH, K)
- Equilibrium composition calculation using Gibbs energy minimization
- Temperature-dependent equilibrium simulations
- Reaction kinetics simulation with Arrhenius rate equations
- Mixture property calculations

## Usage

### Hierarchical Thermodynamic Calculator

```julia
# Run the hierarchical calculator
julia hierarchical_calculator.jl
```

This will:
1. Calculate thermodynamic properties for N2 using different methods
2. Compare theoretical estimations, NASA polynomials, and the hierarchical approach
3. Generate NASA-7 polynomials for both temperature ranges
4. Create visualizations of the calculated properties
5. Demonstrate reaction equilibrium calculations for methane combustion

### Reaction Simulator

```julia
# Run the reaction simulator
julia reaction_simulator.jl
```

This will demonstrate:
1. Parsing chemical reaction equations
2. Calculating reaction thermodynamics (ΔG, ΔH, K)
3. Finding equilibrium compositions at different temperatures
4. Simulating reaction kinetics
5. Generating plots of the results

## API Reference

### Hierarchical Calculator

- `calculate_theoretical_properties(species_name, formula, temperature)`: Calculate properties using theoretical methods
- `progressively_refine_thermodynamic_data(species_name, formula, temperature, data_sources)`: Refine thermodynamic data through data source hierarchy
- `calculate_properties_range(species_name, formula, temp_range, step, data_sources)`: Calculate properties over a temperature range
- `fit_nasa7_coefficients(temps, cp_vals, h_vals, s_vals)`: Fit NASA-7 polynomials to thermodynamic data

### Reaction Simulator

- `parse_reaction_equation(equation)`: Parse a chemical reaction equation
- `calculate_delta_g_reaction(reactants, products, temperature, data_sources)`: Calculate Gibbs energy change for a reaction
- `calculate_delta_h_reaction(reactants, products, temperature, data_sources)`: Calculate enthalpy change for a reaction
- `calculate_equilibrium_constant(delta_g, temperature)`: Calculate equilibrium constant from Gibbs energy
- `calculate_reaction_equilibrium(equation, temperature, pressure, initial_moles, data_sources)`: Calculate equilibrium composition
- `simulate_batch_reaction(equation, initial_composition, time_span, temperature, a_forward, e_act_forward)`: Simulate reaction kinetics
- `simulate_temperature_dependence(equation, initial_composition, temp_range, temp_step, pressure, data_sources)`: Simulate temperature effects on equilibrium

## Data Sources

The hierarchical approach uses a priority-based system for data sources:

1. Theoretical estimates (lowest priority)
2. GRI-MECH
3. CHEMKIN
4. NASA CEA
5. JANAF
6. ThermoML
7. TDE
8. Burcat
9. ATcT (highest priority)

Each source refines previous calculations with proper uncertainty handling.

## Examples

Example output for N2 at 1000K:

```
Cp = 29.58 ± 0.59 J/mol/K
H = 23.45 ± 0.47 kJ/mol
S = 208.31 ± 4.17 J/mol/K
G = -184.86 ± 3.70 kJ/mol
```

NASA-7 polynomial for N2:

```
NASA-7 Polynomial Coefficients for N2 (N2):
Temperature ranges: 200.00-1000.00 K and 1000.00-2000.00 K

! NASA-7 polynomial in standard format
N2              N2          200.00    2000.00    1000.00
 2.95257637e+00 1.39690040e-03-4.92631603e-07 7.86010195e-11-4.60755204e-15    2
-9.23948688e+02 5.87188762e+00 3.53100528e+00-1.23660988e-04-5.02999433e-07    3
 2.43530612e-09-1.40881235e-12-1.04697628e+03 2.96747038e+00                   4
```

Example reaction equilibrium for methane combustion (CH4 + 2O2 = CO2 + 2H2O) at 1000K:

```
Thermodynamic analysis:
  ΔG° = -801.25 ± 16.03 kJ/mol
  K = 7.7664e+041

Equilibrium composition (moles):
  CH4: 0.000000
  CO2: 1.000000
  H2O: 2.000000
  O2: 0.000000
  Reaction extent: 1.000000
```

## Visualization

The package generates several visualizations:

1. Thermodynamic properties (Cp, H, S, G) with uncertainty bands
2. Comparison of different calculation methods
3. Reaction equilibrium as a function of temperature
4. Reaction kinetics over time

## License

MIT License