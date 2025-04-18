# Hierarchical Thermodynamic Data Sources

This table is organized in ascending order of accuracy and reliability. When searching for thermodynamic data for a specific species, start from the top and proceed downward, with each subsequent source potentially overwriting values from previous sources.

| Rank | Database | Species Count | Polynomial Type | Year | Maintainer | Key Characteristics |
|------|----------|---------------|----------------|------|------------|---------------------|
| 1 | GRI-MECH | Limited (~53) | NASA 7 | 1999 (3.0) | Berkeley/Stanford/SRI | Focused on natural gas combustion; small but well-validated |
| 2 | CHEMKIN | ~150 standard + user additions | NASA 7 | Various versions | Originally Sandia, now ANSYS | Basic library for combustion modeling; widely compatible |
| 3 | NASA CEA | 2000+ | NASA 7 & 9 | Latest update 2002, ongoing | NASA Glenn Research Center | Original developer of NASA polynomials; aerospace focus |
| 4 | JANAF | 1800+ | Tabular data | Fourth Edition: 1998 | NIST | Historical gold standard; basis for many other databases |
| 5 | ThermoML | Growing collection | XML-based | Ongoing | NIST/IUPAC | Standardized format for experimental data with uncertainties |
| 6 | ThermoData Engine (TDE) | Millions of data points | Various | 2009, Ongoing | NIST | Dynamically evaluates data; incorporates prediction methods |
| 7 | Burcat Database | 3000+ | NASA 7 & 9 | Latest update 2023 | ReSpecTh/ELTE | Superior in comparative studies; regularly updated |
| 8 | Active Thermochemical Tables (ATcT) | Constantly growing | Various | Ongoing | Argonne National Laboratory | Highest accuracy; ensures thermodynamic consistency |

## Hierarchical Data Search Strategy

1. **Begin with specialized databases** if your application matches their focus (e.g., GRI-MECH for natural gas combustion).
   
2. **Progress to general-purpose databases** (CHEMKIN, NASA CEA) which offer broader coverage but may have older or less refined data for some species.
   
3. **Check NIST resources** (JANAF, ThermoML, TDE) which offer rigorously evaluated data with uncertainty information.
   
4. **Finish with the most accurate sources** (Burcat, ATcT) which feature the latest updates and highest reliability.

## Data Integration Best Practices

- When a species appears in multiple databases, values from higher-ranked sources should generally supersede those from lower-ranked sources.
  
- For critical applications, cross-validate data between multiple high-ranking sources.
  
- Consider temperature ranges carefully - some databases may excel in specific temperature regimes.
  
- For species used in kinetic mechanisms, maintain thermodynamic consistency by using data from the same source for all species in a reaction when possible.
  
- Preserve enthalpy reference states and standard states when mixing data from different sources.
