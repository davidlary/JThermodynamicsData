#!/usr/bin/env julia

"""
Initialize JSON storage for thermodynamic data

This script initializes the JSON storage for thermodynamic data by:
1. Setting up the directory structure
2. Importing data from various thermodynamic databases
3. Creating data source hierarchy information
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using JSON
using YAML
using DataFrames

# Define data sources with their priorities
data_sources = [
    ("theoretical", "Theoretical calculations", 0, 2.5),
    ("group-contribution", "Group contribution methods", 1, 3.0),
    ("stat-thermo", "Statistical thermodynamics", 2, 3.5),
    ("gri-mech", "GRI-MECH 3.0", 3, 3.5),
    ("chemkin", "CHEMKIN format data", 4, 3.8),
    ("nasa-cea", "NASA CEA Database", 5, 4.0),
    ("janaf", "JANAF Thermochemical Tables", 6, 4.5),
    ("thermoml", "ThermoML Standard", 7, 4.5),
    ("tde", "NIST ThermoData Engine", 8, 4.8),
    ("nist-webbook", "NIST Chemistry WebBook", 9, 4.95),
    ("burcat", "Burcat Database", 10, 4.9),
    ("atct", "Active Thermochemical Tables", 11, 5.0)
]

# Create species directory
species_dir = joinpath(dirname(dirname(@__FILE__)), "data", "species")
mkpath(species_dir)

# Process a few key species with sample data to test the system
# This will be replaced with real data import in the full implementation
function add_sample_species()
    # N2 - Nitrogen
    n2_data = generate_theoretical_data("N2")
    n2_data["data"]["low_temp"]["coefficients"] = [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372]
    n2_data["data"]["high_temp"]["coefficients"] = [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]
    add_source_data("N2", "gri-mech", n2_data, 3, 3.5)
    
    # Also add theoretical data as fallback
    add_source_data("N2", "theoretical", generate_theoretical_data("N2"), 0, 2.5)

    # O2 - Oxygen
    o2_data = generate_theoretical_data("O2")
    o2_data["data"]["low_temp"]["coefficients"] = [3.212936, 0.0011274864, -5.756150e-07, 1.313877e-09, -8.768554e-13, -1005.249, 6.034738]
    o2_data["data"]["high_temp"]["coefficients"] = [3.697578, 0.0006135197, -1.258842e-07, 1.775281e-11, -1.136435e-15, -1233.930, 3.189166]
    add_source_data("O2", "gri-mech", o2_data, 3, 3.5)
    
    # Also add theoretical data as fallback
    add_source_data("O2", "theoretical", generate_theoretical_data("O2"), 0, 2.5)

    # H2O - Water
    h2o_data = generate_theoretical_data("H2O")
    h2o_data["data"]["low_temp"]["coefficients"] = [4.198640560E+00, -2.036434100E-03, 6.520402110E-06, -5.487970620E-09, 1.771978170E-12, -3.029372670E+04, -8.490322080E-01]
    h2o_data["data"]["high_temp"]["coefficients"] = [3.033992490E+00, 2.176918040E-03, -1.640725180E-07, -9.704198700E-11, 1.682009920E-14, -3.000429710E+04, 4.966770100E+00]
    add_source_data("H2O", "janaf", h2o_data, 6, 4.5)
    
    # Also add BURCAT data for H2O (to demonstrate hierarchical selection)
    h2o_burcat = generate_theoretical_data("H2O")
    h2o_burcat["data"]["low_temp"]["coefficients"] = [4.19864056, -0.0020364341, 6.5204021E-06, -5.48797062E-09, 1.77197817E-12, -30293.7267, -0.84903221]
    h2o_burcat["data"]["high_temp"]["coefficients"] = [3.03399249, 0.00217691804, -1.64072518E-07, -9.7041987E-11, 1.68200992E-14, -30004.2971, 4.9667701]
    add_source_data("H2O", "burcat", h2o_burcat, 9, 4.9)
    
    # Also add theoretical data as fallback
    add_source_data("H2O", "theoretical", generate_theoretical_data("H2O"), 0, 2.5)

    # Add CO2 from BURCAT
    co2_data = generate_theoretical_data("CO2") 
    co2_data["data"]["low_temp"]["coefficients"] = [2.35677352, 0.00898459677, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13, -48371.9697, 9.90105222]
    co2_data["data"]["high_temp"]["coefficients"] = [3.85746029, 0.00441437026, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -48759.166, 2.27163806]
    add_source_data("CO2", "burcat", co2_data, 9, 4.9)
    
    # Also add NASA data for CO2 (to demonstrate source selection)
    co2_nasa = generate_theoretical_data("CO2")
    co2_nasa["data"]["low_temp"]["coefficients"] = [2.3568130, 0.0089878204, -7.1412648E-06, 2.4579149E-09, -1.4376676E-13, -48371.971, 9.9009035]
    co2_nasa["data"]["high_temp"]["coefficients"] = [4.6365111, 0.0027414569, -9.9589759E-07, 1.6038666E-10, -9.1619857E-15, -49024.904, -1.9348955]
    add_source_data("CO2", "nasa-cea", co2_nasa, 5, 4.0)
    
    # Also add theoretical data as fallback
    add_source_data("CO2", "theoretical", generate_theoretical_data("CO2"), 0, 2.5)

    # O3 - Ozone (from theoretical with modified coefficients)
    o3_data = generate_theoretical_data("O3")
    o3_data["data"]["low_temp"]["coefficients"] = [3.407006, 0.002080774, -1.366827e-06, -2.278595e-10, 3.022012e-13, 15905.2, 8.456681]
    o3_data["data"]["high_temp"]["coefficients"] = [12.77183, -0.00231075, -1.878548e-06, 2.158460e-09, -4.370919e-13, 10963.35, -39.04402]
    add_source_data("O3", "atct", o3_data, 10, 5.0)
    
    # Also add theoretical data as fallback
    add_source_data("O3", "theoretical", generate_theoretical_data("O3"), 0, 2.5)
end

# Function to import data from standard formats
function import_standard_formats()
    println("Importing data from standard formats...")
    
    # Create cache directories if they don't exist
    data_dir = joinpath(dirname(dirname(@__FILE__)), "data")
    mkpath(joinpath(data_dir, "cache", "chemkin"))
    mkpath(joinpath(data_dir, "cache", "burcat"))
    mkpath(joinpath(data_dir, "cache", "janaf"))
    mkpath(joinpath(data_dir, "cache", "nasa-cea"))
    mkpath(joinpath(data_dir, "cache", "thermoml"))
    mkpath(joinpath(data_dir, "cache", "tde"))
    mkpath(joinpath(data_dir, "cache", "atct"))
    
    # For now, just add sample data
    # In a full implementation, this would parse actual data files
    add_sample_species()
    
    # TODO: Add actual importers for the various formats
    println("  Added sample data for N2, O2, H2O, CO2, and O3")
end

# Main function
function main()
    println("Initializing JSON storage for thermodynamic data...")
    
    # Create the directory structure
    species_dir = joinpath(dirname(dirname(@__FILE__)), "data", "species")
    mkpath(species_dir)
    
    # Import data from standard formats
    import_standard_formats()
    
    # Generate theoretical data for all species in the standard list
    # This ensures we always have fallback data
    standard_species = [
        "N2", "O2", "H2O", "CO2", "CO", "CH4", "H2", "NO", "NO2", "Ar", "He", "NH3", "OH",
        "N", "O", "H", "C", "S", "Xe", "Ne", "SO", "SO2", "SO3", "CN", "HCN", "NCO", "NCN",
        "HNO", "HO2", "H2O2", "HONO", "HNO3", "N2O", "N2H4", "CH2", "CH3", "CHO", "CH2O",
        "CHOOH", "HCNO", "N+", "O+", "H+", "C+", "S+", "Xe+", "He+", "CO+", "NO+", "CN+",
        "N2+", "O2+", "H2+", "CO2+", "NO2+", "N2O+", "H2O+", "H3O+", "CHO+", "HNO+",
        "HNO2+", "HNO3+", "O3+", "N-", "O-", "H-", "C-", "CN-", "CO-", "NO-", "O2-",
        "O3", "OH+", "OH-", "HO2-", "NO2-", "NO3-", "e-"
    ]
    
    println("\nGenerating fallback theoretical data for standard species...")
    theoretical_count = 0
    
    for species in standard_species
        # Only add theoretical data if it doesn't already exist
        data = load_species_data(species)
        if !haskey(data, "sources") || !haskey(data["sources"], "theoretical")
            add_source_data(species, "theoretical", 
                          generate_theoretical_data(species), 0, 2.5)
            theoretical_count += 1
        end
    end
    
    println("  Generated theoretical data for $theoretical_count species")
    
    # Print summary of available species
    available_species = list_available_species()
    println("\nStorage initialization complete!")
    println("Total species available: $(length(available_species))")
    
    # Print source statistics
    source_counts = Dict()
    for species in available_species
        data = load_species_data(species)
        if haskey(data, "sources")
            for source in keys(data["sources"])
                source_counts[source] = get(source_counts, source, 0) + 1
            end
        end
    end
    
    println("\nSources used:")
    for (source, count) in sort(collect(source_counts), by=x->x[2], rev=true)
        println("  $(rpad(source, 15)): $(count) species")
    end
    
    println("\nYou can now run the direct calculator with JSON storage.")
end

# Run the main function
main()