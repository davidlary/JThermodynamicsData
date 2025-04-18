"""
Physical and chemical constants used in thermodynamic calculations.
"""

# Physical constants
const R = 8.31446261815324  # Universal gas constant, J/(mol·K)
const NA = 6.02214076e23    # Avogadro's number, 1/mol
const KB = 1.380649e-23     # Boltzmann constant, J/K
const H = 6.62607015e-34    # Planck's constant, J·s
const C = 299792458.0       # Speed of light, m/s
const ELECTRON_MASS = 9.1093837015e-31  # Electron mass, kg
const AMU = 1.66053906660e-27  # Atomic mass unit, kg

# Standard reference conditions
const P_STANDARD = 101325.0  # Standard pressure, Pa
const T_STANDARD = 298.15    # Standard temperature, K

# Conversion factors
const HARTREE_TO_J = 4.3597447222071e-18  # Hartree to Joules
const CAL_TO_J = 4.184                    # Calorie to Joules
const BTU_TO_J = 1055.06                  # BTU to Joules
const EV_TO_J = 1.602176634e-19           # Electron volt to Joules

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
    "Pb" => 207.2
)

# Standard formation enthalpy of common molecules at 298.15 K (in kJ/mol)
const STD_FORMATION_ENTHALPY = Dict(
    "H2" => 0.0,
    "O2" => 0.0,
    "N2" => 0.0,
    "H2O" => -285.83,
    "CO2" => -393.52,
    "CO" => -110.53,
    "CH4" => -74.87,
    "C2H6" => -84.68,
    "C3H8" => -103.85,
    "C2H4" => 52.47,
    "C2H2" => 226.73,
    "NH3" => -45.90,
    "NO" => 90.29,
    "NO2" => 33.10,
    "SO2" => -296.83,
    "H2S" => -20.63
)

# Elements and their valences for group contribution methods
const ELEMENT_VALENCES = Dict(
    "H" => 1,
    "C" => 4,
    "N" => 3,
    "O" => 2,
    "F" => 1,
    "Si" => 4,
    "P" => 3,
    "S" => 2,
    "Cl" => 1,
    "Br" => 1,
    "I" => 1
)

# Database column name constants
const COL_ID = "id"
const COL_SPECIES_NAME = "species_name"
const COL_FORMULA = "formula"
const COL_CAS_NUMBER = "cas_number"
const COL_DATA_SOURCE = "data_source"
const COL_POLYNOMIAL_TYPE = "polynomial_type"
const COL_TEMP_MIN = "temperature_min"
const COL_TEMP_MAX = "temperature_max"
const COL_RELIABILITY_SCORE = "reliability_score"
const COL_DATA_JSON = "data_json"
const COL_METADATA_JSON = "metadata_json"
const COL_DATE_ADDED = "date_added"
const COL_DATE_MODIFIED = "date_modified"