#!/usr/bin/env julia

"""
Import Data from All Sources

This script imports thermodynamic data from all available sources:
- ATcT (Active Thermochemical Tables) - highest priority
- Burcat Database
- TDE (NIST ThermoData Engine)
- ThermoML Standard
- JANAF Thermochemical Tables
- NASA CEA Database
- CHEMKIN format data
- GRI-MECH 3.0

Data is imported into the JSON storage system with proper source priorities.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using JSON
using YAML
using DataFrames
using DuckDB
using Logging

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

"""
    nasa7_to_json_data(a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low,
                     a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high,
                     t_low, t_mid, t_high)

Convert NASA-7 polynomial coefficients to the JSON format used in our storage system.
"""
function nasa7_to_json_data(a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low,
                         a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high,
                         t_low, t_mid, t_high, uncertainty=0.05)
    return Dict(
        "polynomial_type" => "nasa7",
        "temperature_min" => t_low,
        "temperature_max" => t_high,
        "data" => Dict(
            "low_temp" => Dict(
                "range" => [t_low, t_mid],
                "coefficients" => [a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low]
            ),
            "high_temp" => Dict(
                "range" => [t_mid, t_high],
                "coefficients" => [a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high]
            )
        ),
        "uncertainty" => uncertainty
    )
end

"""
    import_from_atct()

Import data from Active Thermochemical Tables (ATcT).
ATcT has the highest priority.
"""
function import_from_atct()
    println("Importing data from Active Thermochemical Tables (ATcT)...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    atct_dir = joinpath(project_dir, "data", "cache", "atct")
    
    try
        # Check for a sample data file first
        sample_file = joinpath(project_dir, "data", "test_data", "sample_atct.dat")
        
        if isfile(sample_file)
            println("  Using sample ATcT data file")
            species_count = 0
            
            # Parse the file
            atct_data = parse_atct_file(sample_file)
            
            # Add each species to our JSON storage
            for (species_name, species_data) in atct_data
                # Check if we have all required data
                if haskey(species_data, "formula") && 
                   haskey(species_data, "h298") && 
                   haskey(species_data, "s298") && 
                   haskey(species_data, "cp")
                    
                    # Convert to NASA-7 polynomial format
                    # (For real implementation, we'd have a proper function to do this)
                    # Here we're generating a reasonable approximation for the demo
                    
                    # Use h298 and s298 to inform our polynomial coefficients
                    h298 = get(species_data, "h298", 0.0)
                    s298 = get(species_data, "s298", 0.0)
                    cp_values = get(species_data, "cp", Dict())
                    
                    # Create basic coefficients with the appropriate h298 and s298 values
                    a1_low = 3.5
                    a6_low = h298 * 4184 / 8.31446  # Convert kcal/mol to J/mol, then to dimensionless
                    a7_low = s298 / 8.31446 - log(298.15)  # Convert to dimensionless
                    
                    # If we have Cp values, adjust a1 (constant term) to match
                    if haskey(cp_values, "298.15")
                        a1_low = cp_values["298.15"] / 8.31446  # Cp at 298K in R units
                    end
                    
                    # High temperature coefficients
                    a1_high = a1_low * 1.05  # Slight increase for high T
                    a6_high = a6_low
                    a7_high = a7_low
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a1_low, 0.0, 0.0, 0.0, 0.0, a6_low, a7_low,
                        a1_high, 0.0, 0.0, 0.0, 0.0, a6_high, a7_high,
                        200.0, 1000.0, 6000.0,
                        0.03  # ATcT typically has low uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "atct", nasa7_data, 10, 5.0)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from ATcT")
        else
            println("  No ATcT data file found. To add real ATcT data, place a data file in data/cache/atct/")
        end
    catch e
        @error "Error importing ATcT data: $e"
    end
end

"""
    import_from_burcat()

Import data from the Burcat thermodynamic database.
Burcat is the second highest priority source.
"""
function import_from_burcat()
    println("Importing data from Burcat Database...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    burcat_dir = joinpath(project_dir, "data", "cache", "burcat")
    burcat_file = joinpath(burcat_dir, "burcat", "BURCAT.THR")
    
    # Check for test data if real data not available
    sample_file = joinpath(project_dir, "data", "test_data", "sample_burcat.dat")
    
    try
        # Determine which file to use
        file_to_use = ""
        if isfile(burcat_file)
            file_to_use = burcat_file
            println("  Using Burcat database file")
        elseif isfile(sample_file)
            file_to_use = sample_file
            println("  Using sample Burcat data file")
        else
            println("  No Burcat data file found. To add real Burcat data, place BURCAT.THR in data/cache/burcat/burcat/")
            return
        end
        
        # Parse the file
        burcat_data = parse_burcat_file(file_to_use)
        species_count = 0
        
        # Add each species to our JSON storage
        for (species_name, species_data) in burcat_data
            if haskey(species_data, "formula") && 
               haskey(species_data, "polynomial_data") &&
               haskey(species_data["polynomial_data"], "a_low") &&
               haskey(species_data["polynomial_data"], "a_high")
                
                # Extract NASA-7 coefficients
                a_low = species_data["polynomial_data"]["a_low"]
                a_high = species_data["polynomial_data"]["a_high"]
                t_ranges = species_data["polynomial_data"]["temperature_ranges"]
                
                # Create NASA-7 polynomial JSON structure
                nasa7_data = nasa7_to_json_data(
                    a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], a_low[7],
                    a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], a_high[7],
                    t_ranges[1], t_ranges[2], t_ranges[3],
                    0.05  # Burcat typically has higher uncertainty than ATcT
                )
                
                # Add to our storage
                add_source_data(species_name, "burcat", nasa7_data, 9, 4.9)
                species_count += 1
            end
        end
        
        println("  Added $(species_count) species from Burcat")
    catch e
        @error "Error importing Burcat data: $e"
    end
end

"""
    import_from_janaf()

Import data from JANAF Thermochemical Tables.
"""
function import_from_janaf()
    println("Importing data from JANAF Thermochemical Tables...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    janaf_dir = joinpath(project_dir, "data", "cache", "janaf")
    
    # Check for test data
    sample_file = joinpath(project_dir, "data", "test_data", "sample_janaf.dat")
    
    try
        if isfile(sample_file)
            println("  Using sample JANAF data file")
            
            # Parse the file
            janaf_data = parse_janaf_file(sample_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in janaf_data
                # Check if we can convert to NASA-7 polynomial
                nasa7_data = janaf_to_nasa7(species_data)
                
                if nasa7_data !== nothing
                    # Extract coefficients
                    a_low = nasa7_data["a_low"]
                    a_high = nasa7_data["a_high"]
                    t_ranges = nasa7_data["temperature_ranges"]
                    
                    # Create JSON data
                    json_data = nasa7_to_json_data(
                        a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], a_low[7],
                        a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], a_high[7],
                        t_ranges[1], t_ranges[2], t_ranges[3],
                        0.07  # JANAF typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "janaf", json_data, 6, 4.5)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from JANAF")
        else
            println("  No JANAF data file found. To add real JANAF data, place data files in data/cache/janaf/")
        end
    catch e
        @error "Error importing JANAF data: $e"
    end
end

"""
    import_from_nasa_cea()

Import data from NASA CEA Database.
"""
function import_from_nasa_cea()
    println("Importing data from NASA CEA Database...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    nasa_dir = joinpath(project_dir, "data", "cache", "nasa-cea")
    nasa_file = joinpath(nasa_dir, "thermo.inp")
    
    # Check for a test file
    sample_file = joinpath(project_dir, "data", "test_data", "sample_nasa7.dat")
    
    try
        # Determine which file to use
        file_to_use = ""
        if isfile(nasa_file)
            file_to_use = nasa_file
            println("  Using NASA CEA database file")
        elseif isfile(sample_file)
            file_to_use = sample_file
            println("  Using sample NASA-7 data file")
        else
            println("  No NASA CEA data file found. To add real NASA data, place thermo.inp in data/cache/nasa-cea/")
            return
        end
        
        # Parse the file (NASA format is similar to CHEMKIN)
        nasa_data = parse_nasa7_file(file_to_use)
        species_count = 0
        
        # Add each species to our JSON storage
        for (species_name, species_data) in nasa_data
            if haskey(species_data, "polynomial_data") &&
               haskey(species_data["polynomial_data"], "a_low") &&
               haskey(species_data["polynomial_data"], "a_high")
                
                # Extract NASA-7 coefficients
                a_low = species_data["polynomial_data"]["a_low"]
                a_high = species_data["polynomial_data"]["a_high"]
                t_ranges = species_data["polynomial_data"]["temperature_ranges"]
                
                # Create NASA-7 polynomial JSON structure
                nasa7_data = nasa7_to_json_data(
                    a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], a_low[7],
                    a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], a_high[7],
                    t_ranges[1], t_ranges[2], t_ranges[3],
                    0.08  # NASA-CEA typical uncertainty
                )
                
                # Add to our storage
                add_source_data(species_name, "nasa-cea", nasa7_data, 5, 4.0)
                species_count += 1
            end
        end
        
        println("  Added $(species_count) species from NASA CEA")
    catch e
        @error "Error importing NASA CEA data: $e"
    end
end

"""
    import_from_chemkin()

Import data from CHEMKIN format files.
"""
function import_from_chemkin()
    println("Importing data from CHEMKIN format files...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    chemkin_dir = joinpath(project_dir, "data", "cache", "chemkin")
    
    # Check for test data
    sample_file = joinpath(project_dir, "data", "test_data", "sample_chemkin.dat")
    
    try
        if isfile(sample_file)
            println("  Using sample CHEMKIN data file")
            
            # Parse the file
            chemkin_data = parse_chemkin_file(sample_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in chemkin_data
                if haskey(species_data, "nasa7") &&
                   haskey(species_data["nasa7"], "a_low") &&
                   haskey(species_data["nasa7"], "a_high")
                    
                    # Extract NASA-7 coefficients
                    a_low = species_data["nasa7"]["a_low"]
                    a_high = species_data["nasa7"]["a_high"]
                    t_ranges = species_data["nasa7"]["temperature_ranges"]
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], a_low[7],
                        a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], a_high[7],
                        t_ranges[1], t_ranges[2], t_ranges[3],
                        0.10  # CHEMKIN typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "chemkin", nasa7_data, 4, 3.8)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from CHEMKIN")
        else
            println("  No CHEMKIN data file found. To add real CHEMKIN data, place files in data/cache/chemkin/")
        end
    catch e
        @error "Error importing CHEMKIN data: $e"
    end
end

"""
    import_from_gri_mech()

Import data from GRI-MECH 3.0.
GRI-MECH is a specific CHEMKIN-format mechanism for natural gas combustion.
"""
function import_from_gri_mech()
    println("Importing data from GRI-MECH 3.0...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    gri_file = joinpath(project_dir, "data", "cache", "gri-mech", "therm.dat")
    
    try
        if isfile(gri_file)
            println("  Using GRI-MECH 3.0 data file")
            
            # Parse the file (GRI-MECH is in CHEMKIN format)
            gri_data = parse_chemkin_file(gri_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in gri_data
                if haskey(species_data, "nasa7") &&
                   haskey(species_data["nasa7"], "a_low") &&
                   haskey(species_data["nasa7"], "a_high")
                    
                    # Extract NASA-7 coefficients
                    a_low = species_data["nasa7"]["a_low"]
                    a_high = species_data["nasa7"]["a_high"]
                    t_ranges = species_data["nasa7"]["temperature_ranges"]
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6], a_low[7],
                        a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6], a_high[7],
                        t_ranges[1], t_ranges[2], t_ranges[3],
                        0.10  # GRI-MECH typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "gri-mech", nasa7_data, 3, 3.5)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from GRI-MECH")
        else
            println("  No GRI-MECH data file found. To add GRI-MECH data, download therm.dat to data/cache/gri-mech/")
        end
    catch e
        @error "Error importing GRI-MECH data: $e"
    end
end

"""
    import_from_nist_webbook()

Import data from NIST Chemistry WebBook.
NIST WebBook is a high priority source with reliable reference data.
"""
function import_from_nist_webbook()
    println("Importing data from NIST Chemistry WebBook...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    nist_dir = joinpath(project_dir, "data", "cache", "nist-webbook")
    mkpath(nist_dir)
    
    try
        # Check for test data
        sample_file = joinpath(project_dir, "data", "test_data", "sample_nist.json")
        
        if isfile(sample_file)
            println("  Using sample NIST WebBook data file")
            
            # Parse the file (assuming JSON format for simplicity)
            nist_data = JSON.parsefile(sample_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in nist_data
                if haskey(species_data, "properties") &&
                   haskey(species_data["properties"], "heat_capacity")
                    
                    # Get a representative heat capacity value
                    cp_data = species_data["properties"]["heat_capacity"]
                    cp_298 = get(cp_data, "298.15", 30.0) / 8.31446  # Convert to R units
                    
                    # Get enthalpy and entropy at 298K if available
                    h_298 = 0.0
                    s_298 = 0.0
                    
                    if haskey(species_data["properties"], "enthalpy")
                        h_298 = get(species_data["properties"]["enthalpy"], "298.15", 0.0) / 8.31446
                    end
                    
                    if haskey(species_data["properties"], "entropy")
                        s_298 = get(species_data["properties"]["entropy"], "298.15", 0.0) / 8.31446
                    end
                    
                    # Create a simple NASA-7 polynomial
                    a1_low = cp_298
                    a6_low = h_298
                    a7_low = s_298 - log(298.15)
                    
                    a1_high = cp_298 * 1.05  # Slight increase for high temperatures
                    a6_high = a6_low
                    a7_high = a7_low
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a1_low, 0.0, 0.0, 0.0, 0.0, a6_low, a7_low,
                        a1_high, 0.0, 0.0, 0.0, 0.0, a6_high, a7_high,
                        200.0, 1000.0, 6000.0,
                        0.04  # NIST WebBook typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "nist-webbook", nasa7_data, 9, 4.95)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from NIST WebBook")
        else
            println("  No NIST WebBook data file found. To add real NIST data, create a file at data/cache/nist-webbook/")
        end
    catch e
        @error "Error importing NIST WebBook data: $e"
    end
end

"""
    import_from_tde()

Import data from NIST ThermoData Engine (TDE).
TDE is a high priority source.
"""
function import_from_tde()
    println("Importing data from NIST ThermoData Engine...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    tde_dir = joinpath(project_dir, "data", "cache", "tde")
    
    # Check for test data
    sample_file = joinpath(project_dir, "data", "test_data", "sample_tde.json")
    
    try
        if isfile(sample_file)
            println("  Using sample TDE data file")
            
            # Parse the file
            tde_data = parse_tde_file(sample_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in tde_data
                # TDE data typically includes tabular Cp, H, S data that needs conversion to NASA-7
                # For this demo, we'll create simplified NASA-7 polynomials
                
                if haskey(species_data, "properties") &&
                   haskey(species_data["properties"], "heat_capacity")
                    
                    # Get a representative heat capacity value
                    cp_data = species_data["properties"]["heat_capacity"]
                    cp_298 = get(cp_data, "298.15", 30.0) / 8.31446  # Convert to R units
                    
                    # Get enthalpy and entropy at 298K if available
                    h_298 = 0.0
                    s_298 = 0.0
                    
                    if haskey(species_data["properties"], "enthalpy")
                        h_298 = get(species_data["properties"]["enthalpy"], "298.15", 0.0) / 8.31446
                    end
                    
                    if haskey(species_data["properties"], "entropy")
                        s_298 = get(species_data["properties"]["entropy"], "298.15", 0.0) / 8.31446
                    end
                    
                    # Create a simple NASA-7 polynomial
                    a1_low = cp_298
                    a6_low = h_298
                    a7_low = s_298 - log(298.15)
                    
                    a1_high = cp_298 * 1.1  # Slight increase for high temperatures
                    a6_high = a6_low
                    a7_high = a7_low
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a1_low, 0.0, 0.0, 0.0, 0.0, a6_low, a7_low,
                        a1_high, 0.0, 0.0, 0.0, 0.0, a6_high, a7_high,
                        200.0, 1000.0, 6000.0,
                        0.05  # TDE typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "tde", nasa7_data, 8, 4.8)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from TDE")
        else
            println("  No TDE data file found. To add real TDE data, export files to data/cache/tde/")
        end
    catch e
        @error "Error importing TDE data: $e"
    end
end

"""
    import_from_thermoml()

Import data from ThermoML standard files.
ThermoML is the fourth highest priority source.
"""
function import_from_thermoml()
    println("Importing data from ThermoML standard...")
    
    # Set up the data directory
    project_dir = dirname(dirname(@__FILE__))
    thermoml_dir = joinpath(project_dir, "data", "cache", "thermoml")
    
    # Check for test data
    sample_file = joinpath(project_dir, "data", "test_data", "sample_thermoml.dat")
    
    try
        if isfile(sample_file)
            println("  Using sample ThermoML data file")
            
            # Parse the file
            thermoml_data = parse_thermoml_file(sample_file)
            species_count = 0
            
            # Add each species to our JSON storage
            for (species_name, species_data) in thermoml_data
                # ThermoML data typically includes property measurements that need conversion to NASA-7
                # For this demo, we'll create simplified NASA-7 polynomials
                
                if haskey(species_data, "properties")
                    # Get heat capacity data if available
                    cp_298 = 30.0 / 8.31446  # Default in R units
                    
                    if haskey(species_data["properties"], "heat_capacity")
                        cp_data = species_data["properties"]["heat_capacity"]
                        cp_298 = get(cp_data, "298.15", cp_298)
                    end
                    
                    # Get enthalpy and entropy at 298K if available
                    h_298 = 0.0
                    s_298 = 60.0 / 8.31446  # Default entropy in R units
                    
                    if haskey(species_data["properties"], "enthalpy")
                        h_298 = get(species_data["properties"]["enthalpy"], "298.15", 0.0) / 8.31446
                    end
                    
                    if haskey(species_data["properties"], "entropy")
                        s_298 = get(species_data["properties"]["entropy"], "298.15", s_298) / 8.31446
                    end
                    
                    # Create a simple NASA-7 polynomial
                    a1_low = cp_298
                    a6_low = h_298
                    a7_low = s_298 - log(298.15)
                    
                    a1_high = cp_298 * 1.1  # Slight increase for high temperatures
                    a6_high = a6_low
                    a7_high = a7_low
                    
                    # Create NASA-7 polynomial JSON structure
                    nasa7_data = nasa7_to_json_data(
                        a1_low, 0.0, 0.0, 0.0, 0.0, a6_low, a7_low,
                        a1_high, 0.0, 0.0, 0.0, 0.0, a6_high, a7_high,
                        200.0, 1000.0, 6000.0,
                        0.07  # ThermoML typical uncertainty
                    )
                    
                    # Add to our storage
                    add_source_data(species_name, "thermoml", nasa7_data, 7, 4.5)
                    species_count += 1
                end
            end
            
            println("  Added $(species_count) species from ThermoML")
        else
            println("  No ThermoML data file found. To add real ThermoML data, place files in data/cache/thermoml/")
        end
    catch e
        @error "Error importing ThermoML data: $e"
    end
end

"""
    generate_theoretical_data_for_standard_species()

Generate theoretical data for all standard species as a fallback.
"""
function generate_theoretical_data_for_standard_species()
    println("\nGenerating fallback theoretical data for standard species...")
    
    # Standard list of chemical species
    standard_species = [
        "N2", "O2", "H2O", "CO2", "CO", "CH4", "H2", "NO", "NO2", "Ar", "He", "NH3", "OH",
        "N", "O", "H", "C", "S", "Xe", "Ne", "SO", "SO2", "SO3", "CN", "HCN", "NCO", "NCN",
        "HNO", "HO2", "H2O2", "HONO", "HNO3", "N2O", "N2H4", "CH2", "CH3", "CHO", "CH2O",
        "CHOOH", "HCNO", "N+", "O+", "H+", "C+", "S+", "Xe+", "He+", "CO+", "NO+", "CN+",
        "N2+", "O2+", "H2+", "CO2+", "NO2+", "N2O+", "H2O+", "H3O+", "CHO+", "HNO+",
        "HNO2+", "HNO3+", "O3+", "N-", "O-", "H-", "C-", "CN-", "CO-", "NO-", "O2-",
        "O3", "OH+", "OH-", "HO2-", "NO2-", "NO3-", "e-"
    ]
    
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
end

"""
    print_data_summary()

Print a summary of the available thermodynamic data.
"""
function print_data_summary()
    # Print summary of available species
    available_species = list_available_species()
    println("\nData import complete!")
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
    for (source_name, _, priority, _) in data_sources
        count = get(source_counts, source_name, 0)
        println("  $(rpad(source_name, 20)) (priority: $(priority)): $(count) species")
    end
    
    # Print species with multiple sources
    species_with_multiple_sources = []
    for species in available_species
        data = load_species_data(species)
        if haskey(data, "sources") && length(keys(data["sources"])) > 1
            push!(species_with_multiple_sources, (species, length(keys(data["sources"]))))
        end
    end
    
    if !isempty(species_with_multiple_sources)
        println("\nSpecies with data from multiple sources:")
        for (species, count) in sort(species_with_multiple_sources, by=x->x[2], rev=true)[1:min(10, length(species_with_multiple_sources))]
            sources = list_available_sources(species)
            source_text = join(["$(src) (priority: $(pri))" for (src, pri, _) in sources[1:min(3, length(sources))]], ", ")
            if length(sources) > 3
                source_text *= ", ..."
            end
            println("  $(rpad(species, 12)): $(count) sources - $(source_text)")
        end
    end
    
    println("\nYou can now run the direct calculator with JSON storage to visualize the data.")
    println("The best source will be selected based on the hierarchical priority.")
end

"""
    main()

Main function to import data from all sources.
"""
function main()
    println("Importing thermodynamic data from all sources...")
    
    # Create necessary directories
    project_dir = dirname(dirname(@__FILE__))
    mkpath(joinpath(project_dir, "data", "species"))
    
    # Import data from each source in order of priority
    import_from_atct()      # Highest priority (11)
    import_from_burcat()    # Priority 10
    import_from_nist_webbook() # Priority 9
    import_from_tde()       # Priority 8
    import_from_thermoml()  # Priority 7
    import_from_janaf()     # Priority 6
    import_from_nasa_cea()  # Priority 5
    import_from_chemkin()   # Priority 4
    import_from_gri_mech()  # Priority 3 (lowest experimental before theoretical)
    
    # Generate theoretical fallback data
    generate_theoretical_data_for_standard_species()
    
    # Print summary
    print_data_summary()
end

# Run the main function
main()