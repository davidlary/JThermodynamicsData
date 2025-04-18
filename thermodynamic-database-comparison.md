# Comparison of Thermodynamic Data Sources

| Database | Species Count | Polynomial Type | Year | Reliability | Maintainer | Special Features |
|----------|---------------|----------------|------|-------------|------------|-----------------|
| NASA CEA | 2000+ | NASA 7 & 9 | Latest update 2002, ongoing | High - Extensively validated in aerospace applications | NASA Glenn Research Center | Includes gaseous and condensed species; original developer of NASA polynomials format |
| Burcat Database | 3000+ | NASA 7 & 9 | Latest update 2023 | Very High - Considered superior in comparative studies | Moved to ReSpecTh/ELTE (originally Technion) | Extensively used for combustion applications; constantly updated; available in multiple formats |
| JANAF | 1800+ | Tabular data | First: 1965, Fourth Edition: 1998 | High - Gold standard for many years | NIST | Original source of many thermodynamic values; available as PDF tables |
| ThermoData Engine (TDE) | Millions of data points | Various | First: 2009, Ongoing | High - NIST curated | NIST | Dynamically generates recommended data; incorporates property prediction schemes |
| Active Thermochemical Tables (ATcT) | Constantly growing | Various | Ongoing | Very High - Uses network analysis for consistency | Argonne National Laboratory | Uses network approach to ensure thermodynamic consistency across multiple sources |
| CHEMKIN | ~150 standard + user additions | NASA 7 | Various versions | Moderate to High | Originally Sandia, now ANSYS | Widely used in combustion modeling; designed for easy integration with reaction mechanisms |
| GRI-MECH | Limited (~53) | NASA 7 | 1999 (3.0) | High for included species | Berkeley/Stanford/SRI | Focused on natural gas combustion; carefully validated subset |
| ThermoML | Growing collection | XML-based | Ongoing | High - Peer-reviewed data | NIST/IUPAC | Machine-readable format; includes experimental uncertainty data |

## Notes on Reliability and Selection

1. **Reliability Assessment**:
   - The Burcat database is often cited as having superior features when compared to other databases 
   - NIST databases (JANAF, TDE) are considered highly reliable due to their rigorous evaluation processes
   - ATcT uses innovative statistical approaches to ensure thermodynamic consistency

2. **Choosing a Database**:
   - For combustion applications: Burcat or NASA CEA are most comprehensive
   - For general reference with uncertainty data: ThermoData Engine
   - For absolute consistency across species: Active Thermochemical Tables
   - For specialized applications: Domain-specific databases like GRI-MECH (natural gas)

3. **Format Considerations**:
   - NASA 7-coefficient polynomials: Most widely supported in simulation tools
   - NASA 9-coefficient polynomials: Better accuracy, especially at temperature extremes
   - ThermoML: Best for data exchange and preservation of experimental details

4. **Database Integration**:
   - Most modern simulation software can import NASA polynomial formats
   - Many databases provide conversion utilities between formats
   - Most comprehensive approach is to use multiple sources with cross-validation
