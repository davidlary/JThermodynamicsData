# Thermodynamic Data Sources for O2

## Species Information
- **Species Name**: O2
- **Formula**: O2
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
| Cp | 34.8247 | ±0.7114 | J/mol/K |
| H | 428.8504 | ±1282.5947 | kJ/mol |
| S | 244.4939 | ±5.4667 | J/mol/K |
| G | 184.3565 | ±1280.5007 | kJ/mol |

## Metadata
- **Generated**: 2025-04-18 13:04:33
- **JThermodynamicsData Version**: 1.0.0
