"""
Core type definitions for the JThermodynamicsData package.
"""

export ThermodynamicProperty, PolynomialType, DataSource
export ThermodynamicData, UncertaintyMethod, DatabaseEntry

"""
    PolynomialType

Enum for different polynomial types used in thermodynamic data.
"""
@enum PolynomialType begin
    NASA7
    NASA9
    SHOMATE
    WILHOIT
    TABULAR
    CUSTOM
end

"""
    DataSource

Enum for different data sources.
"""
@enum DataSource begin
    GRIMECH
    CHEMKIN
    NASA_CEA
    JANAF
    THERMOML
    THERMODATA_ENGINE
    BURCAT
    ATCT
    THEORETICAL_GROUP_CONTRIBUTION
    THEORETICAL_STATISTICAL_THERMO
    THEORETICAL_MACHINE_LEARNING
    THEORETICAL_QUANTUM_CHEMISTRY
    UNKNOWN
end

"""
    UncertaintyMethod

Enum for different uncertainty quantification methods.
"""
@enum UncertaintyMethod begin
    MEASUREMENTS
    MONTE_CARLO
    INTERVAL
    NONE
end

"""
    ThermodynamicProperty

Abstract type for thermodynamic properties.
"""
abstract type ThermodynamicProperty end

"""
    CpProperty <: ThermodynamicProperty

Heat capacity at constant pressure.
"""
struct CpProperty <: ThermodynamicProperty
    temperature::Union{Float64, Measurement{Float64}, Particles{Float64}}
    value::Union{Float64, Measurement{Float64}, Particles{Float64}}
    units::String
end

"""
    EnthalpyProperty <: ThermodynamicProperty

Enthalpy property.
"""
struct EnthalpyProperty <: ThermodynamicProperty
    temperature::Union{Float64, Measurement{Float64}, Particles{Float64}}
    value::Union{Float64, Measurement{Float64}, Particles{Float64}}
    units::String
end

"""
    EntropyProperty <: ThermodynamicProperty

Entropy property.
"""
struct EntropyProperty <: ThermodynamicProperty
    temperature::Union{Float64, Measurement{Float64}, Particles{Float64}}
    value::Union{Float64, Measurement{Float64}, Particles{Float64}}
    units::String
end

"""
    GibbsProperty <: ThermodynamicProperty

Gibbs free energy property.
"""
struct GibbsProperty <: ThermodynamicProperty
    temperature::Union{Float64, Measurement{Float64}, Particles{Float64}}
    value::Union{Float64, Measurement{Float64}, Particles{Float64}}
    units::String
end

"""
    PolynomialCoefficients

Structure to hold polynomial coefficients used for calculating thermodynamic properties.
"""
struct PolynomialCoefficients
    type::PolynomialType
    temperature_ranges::Vector{Tuple{Float64, Float64}}
    coefficients::Vector{Vector{Union{Float64, Measurement{Float64}, Particles{Float64}}}}
    source::DataSource
    reliability_score::Float64
    date_modified::String
end

"""
    ThermodynamicData

Structure to hold all thermodynamic data for a species.
"""
mutable struct ThermodynamicData
    species_name::String
    formula::String
    cas_number::String
    molecular_weight::Union{Float64, Measurement{Float64}, Particles{Float64}}
    polynomials::Vector{PolynomialCoefficients}
    tabular_data::Union{Dict{String, Any}, Nothing}
    source_priority::Int
    uncertainty_method::UncertaintyMethod
    properties_cache::Dict{String, Any}
    metadata::Dict{String, Any}
    
    # Constructor with sensible defaults
    function ThermodynamicData(
        species_name::String,
        formula::String,
        cas_number::String = "",
        molecular_weight = 0.0,
        polynomials = PolynomialCoefficients[],
        tabular_data = nothing,
        source_priority = 0,
        uncertainty_method = NONE,
        properties_cache = Dict{String, Any}(),
        metadata = Dict{String, Any}()
    )
        return new(
            species_name, formula, cas_number, molecular_weight,
            polynomials, tabular_data, source_priority, uncertainty_method,
            properties_cache, metadata
        )
    end
end

"""
    DatabaseEntry

Structure representing a row in the thermodynamics database.
"""
struct DatabaseEntry
    id::Int
    species_name::String
    formula::String
    cas_number::String
    data_source::DataSource
    polynomial_type::PolynomialType
    temperature_min::Float64
    temperature_max::Float64
    reliability_score::Float64
    data_json::String
    metadata_json::String
    date_added::String
    date_modified::String
end

"""
    load_data_source(source_str::String)
    
Convert string to DataSource enum.
"""
function load_data_source(source_str::String)
    source_map = Dict(
        "gri-mech" => GRIMECH,
        "chemkin" => CHEMKIN,
        "nasa-cea" => NASA_CEA,
        "janaf" => JANAF,
        "thermoml" => THERMOML,
        "thermodata-engine" => THERMODATA_ENGINE,
        "burcat" => BURCAT,
        "atct" => ATCT,
        "group_contribution" => THEORETICAL_GROUP_CONTRIBUTION,
        "statistical_thermodynamics" => THEORETICAL_STATISTICAL_THERMO,
        "machine_learning" => THEORETICAL_MACHINE_LEARNING,
        "quantum_chemistry" => THEORETICAL_QUANTUM_CHEMISTRY
    )
    
    return get(source_map, lowercase(source_str), UNKNOWN)
end

"""
    load_polynomial_type(type_str::String)
    
Convert string to PolynomialType enum.
"""
function load_polynomial_type(type_str::String)
    type_map = Dict(
        "nasa7" => NASA7,
        "nasa9" => NASA9,
        "shomate" => SHOMATE,
        "wilhoit" => WILHOIT,
        "tabular" => TABULAR,
        "custom" => CUSTOM
    )
    
    return get(type_map, lowercase(type_str), CUSTOM)
end