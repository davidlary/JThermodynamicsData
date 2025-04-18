# Thermodynamic Data Sources for N2

## Species Information
- **Species Name**: N2
- **Formula**: N2
- **Temperature Range**: 100-10000 K

## Data Sources Used
The following data sources were used in hierarchical refinement, listed in order of application:

| Priority | Source | Reliability Score | Weight Factor |
|----------|--------|-------------------|---------------|
| 0 | THEORETICAL | 2.5 | 1.0 |
| 1 | GRI-MECH | 3.5 | 0.517 |
| 3 | NASA-CEA | 4.0 | 0.624 |
| 7 | Burcat | 4.7 | 0.8398 |

## Hierarchical Refinement Process
The thermodynamic properties for this species were calculated using a hierarchical approach:

1. Started with theoretical estimates (statistical thermodynamics and group contribution methods)
2. Progressively refined with experimental data sources in order of priority
3. Applied weighted averaging based on reliability scores and priority levels
4. Propagated uncertainties throughout the refinement process

## Representative Values at 1000.0 K
| Property | Value | Uncertainty | Units |
|----------|-------|-------------|-------|
| Cp | 32.6773 | ±0.6481 | J/mol/K |
| H | 427.6405 | ±1282.6791 | kJ/mol |
| S | 229.3291 | ±5.5783 | J/mol/K |
| G | 198.3114 | ±1280.0567 | kJ/mol |

## Metadata
- **Generated**: 2025-04-18 13:04:33
- **JThermodynamicsData Version**: 1.0.0
