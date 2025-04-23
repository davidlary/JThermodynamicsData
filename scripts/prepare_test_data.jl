#!/usr/bin/env julia

"""
Prepare Test Data Files

This script prepares test data files for all thermodynamic sources.
It creates sample data files in the proper locations without trying
to download from external sources, which may be unreliable.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JSON
using YAML
using Dates

"""
    create_cache_directories()

Create all required cache directories for data files.
"""
function create_cache_directories()
    println("Creating cache directories...")
    
    # Root project directory
    project_dir = dirname(dirname(@__FILE__))
    
    # Create cache directory structure
    for source in ["atct", "burcat", "nist-webbook", "tde", "thermoml", 
                   "janaf", "nasa-cea", "chemkin", "gri-mech"]
        cache_dir = joinpath(project_dir, "data", "cache", source)
        mkpath(cache_dir)
        println("  Created directory: $(cache_dir)")
    end
    
    println("Cache directories created.")
end

"""
    prepare_sample_data(source_name, file_name, content)

Prepare a sample data file for a source.
"""
function prepare_sample_data(source_name, file_name, content)
    # Root project directory
    project_dir = dirname(dirname(@__FILE__))
    
    # Create directory
    cache_dir = joinpath(project_dir, "data", "cache", source_name)
    mkpath(cache_dir)
    
    # Write the file
    file_path = joinpath(cache_dir, file_name)
    
    open(file_path, "w") do f
        write(f, content)
    end
    
    println("  Created sample data file: $(file_path)")
    return file_path
end

"""
    copy_sample_file(source_name, target_filename)

Copy a sample file from test_data to the cache directory.
"""
function copy_sample_file(source_name, source_filename, target_filename)
    # Root project directory
    project_dir = dirname(dirname(@__FILE__))
    
    # Source file
    sample_file = joinpath(project_dir, "data", "test_data", source_filename)
    
    if isfile(sample_file)
        # Create target directory
        cache_dir = joinpath(project_dir, "data", "cache", source_name)
        mkpath(cache_dir)
        
        # Copy file
        target_file = joinpath(cache_dir, target_filename)
        cp(sample_file, target_file, force=true)
        
        println("  Copied sample file $(source_filename) to $(target_file)")
        return true
    else
        println("  Warning: Sample file $(sample_file) not found")
        return false
    end
end

"""
    prepare_atct_data()

Prepare ATcT (Active Thermochemical Tables) sample data.
"""
function prepare_atct_data()
    println("Preparing ATcT data...")
    
    # Copy sample file
    if copy_sample_file("atct", "sample_atct.dat", "ATcT_data.txt")
        println("ATcT data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = """
        # Active Thermochemical Tables data (sample)
        # Species | Formula | H298 (kJ/mol) | S298 (J/mol/K) | Cp (J/mol/K)
        O3 | O3 | 142.7 | 238.93 | 39.199
        H2O | H2O | -285.83 | 69.95 | 33.577
        O2 | O2 | 0.0 | 205.08 | 29.355
        N2 | N2 | 0.0 | 191.61 | 29.125
        CO2 | CO2 | -393.51 | 213.79 | 37.12
        """
        
        prepare_sample_data("atct", "ATcT_data.txt", sample_data)
        println("ATcT sample data created.")
        return true
    end
end

"""
    prepare_burcat_data()

Prepare Burcat database sample data.
"""
function prepare_burcat_data()
    println("Preparing Burcat database...")
    
    # Create burcat subdirectory
    project_dir = dirname(dirname(@__FILE__))
    mkpath(joinpath(project_dir, "data", "cache", "burcat", "burcat"))
    
    # Copy sample file
    if copy_sample_file("burcat/burcat", "sample_burcat.dat", "BURCAT.THR")
        println("Burcat data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = """
        THERMO
           300.      1000.      5000.
        N2               121386N   2               G  200.000  6000.000  1000.00      1
         2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
        -9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
         2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00                   4
        O2               121386O   2               G  200.000  6000.000  1000.00      1
         3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
        -1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
        -9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00                   4
        H2O               20387H   2O   1          G  200.000  6000.000  1000.00      1
         2.67703787E+00 2.97318329E-03-7.73769690E-07 9.44336689E-11-4.26900959E-15    2
        -2.98858938E+04 6.88255571E+00 4.19863520E+00-2.03640170E-03 6.52034160E-06    3
        -5.48792690E-09 1.77196800E-12-3.02937260E+04-8.49009010E-01                   4
        CO2               121386C   1O   2          G  200.000  6000.000  1000.00      1
         4.63651110E+00 2.74145690E-03-9.95897590E-07 1.60386660E-10-9.16198570E-15    2
        -4.90249040E+04-1.93489550E+00 2.35677520E+00 8.98459680E-03-7.12356270E-06    3
         2.45919220E-09-1.43699550E-13-4.83719710E+04 9.90105220E+00                   4
        O3                121386O   3               G  200.000  6000.000  1000.00      1
         1.23302914E+01-1.19324783E-02 7.98741278E-06-1.77194552E-09 1.26075824E-13    2
         1.26755831E+04-4.08823374E+01 3.40738221E+00 2.05379063E-03-1.38283480E-06    3
        -6.27993680E-10 6.54353945E-13 1.58644979E+04 8.28152238E+00                   4
        END
        """
        
        prepare_sample_data("burcat/burcat", "BURCAT.THR", sample_data)
        println("Burcat sample data created.")
        return true
    end
end

"""
    prepare_nist_webbook_data()

Prepare NIST Chemistry WebBook data.
"""
function prepare_nist_webbook_data()
    println("Preparing NIST Chemistry WebBook data...")
    
    # Copy sample file
    if copy_sample_file("nist-webbook", "sample_nist.json", "nist_webbook_data.json")
        println("NIST WebBook data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = Dict(
            "N2" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 29.124,
                        "500.0" => 29.580,
                        "1000.0" => 32.698
                    ),
                    "enthalpy" => Dict(
                        "298.15" => 0.0,
                        "500.0" => 5.829,
                        "1000.0" => 20.653
                    ),
                    "entropy" => Dict(
                        "298.15" => 191.61,
                        "500.0" => 209.41,
                        "1000.0" => 231.96
                    )
                ),
                "formula" => "N2",
                "molecular_weight" => 28.0134
            ),
            "O2" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 29.376,
                        "500.0" => 31.091,
                        "1000.0" => 34.870
                    ),
                    "enthalpy" => Dict(
                        "298.15" => 0.0,
                        "500.0" => 6.085,
                        "1000.0" => 22.542
                    ),
                    "entropy" => Dict(
                        "298.15" => 205.15,
                        "500.0" => 222.52,
                        "1000.0" => 246.13
                    )
                ),
                "formula" => "O2",
                "molecular_weight" => 31.9988
            )
        )
        
        # Add other key species
        sample_data["H2O"] = Dict(
            "properties" => Dict(
                "heat_capacity" => Dict(
                    "298.15" => 33.590,
                    "500.0" => 35.226,
                    "1000.0" => 41.217
                ),
                "enthalpy" => Dict(
                    "298.15" => -241.826,
                    "500.0" => -238.918,
                    "1000.0" => -226.449
                ),
                "entropy" => Dict(
                    "298.15" => 188.84,
                    "500.0" => 207.42,
                    "1000.0" => 234.90
                )
            ),
            "formula" => "H2O",
            "molecular_weight" => 18.0153
        )
        
        sample_data["CO2"] = Dict(
            "properties" => Dict(
                "heat_capacity" => Dict(
                    "298.15" => 37.129,
                    "500.0" => 44.626,
                    "1000.0" => 54.308
                ),
                "enthalpy" => Dict(
                    "298.15" => -393.522,
                    "500.0" => -389.520,
                    "1000.0" => -371.059
                ),
                "entropy" => Dict(
                    "298.15" => 213.79,
                    "500.0" => 235.46,
                    "1000.0" => 269.30
                )
            ),
            "formula" => "CO2",
            "molecular_weight" => 44.0095
        )
        
        sample_data["O3"] = Dict(
            "properties" => Dict(
                "heat_capacity" => Dict(
                    "298.15" => 39.199,
                    "500.0" => 45.037,
                    "1000.0" => 53.35
                ),
                "enthalpy" => Dict(
                    "298.15" => 142.7,
                    "500.0" => 151.5,
                    "1000.0" => 179.7
                ),
                "entropy" => Dict(
                    "298.15" => 238.93,
                    "500.0" => 260.57,
                    "1000.0" => 292.48
                )
            ),
            "formula" => "O3",
            "molecular_weight" => 47.9982
        )
        
        # Save as JSON
        json_str = JSON.json(sample_data, 4)
        prepare_sample_data("nist-webbook", "nist_webbook_data.json", json_str)
        println("NIST WebBook sample data created.")
        return true
    end
end

"""
    prepare_tde_data()

Prepare TDE (NIST ThermoData Engine) sample data.
"""
function prepare_tde_data()
    println("Preparing NIST ThermoData Engine data...")
    
    # Copy sample file
    if copy_sample_file("tde", "sample_tde.json", "tde_data.json")
        println("TDE data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = Dict(
            "N2" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 29.12,
                        "500.0" => 29.58,
                        "1000.0" => 32.70
                    )
                ),
                "formula" => "N2"
            ),
            "O2" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 29.38,
                        "500.0" => 31.09,
                        "1000.0" => 34.87
                    )
                ),
                "formula" => "O2"
            ),
            "H2O" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 33.59,
                        "500.0" => 35.23,
                        "1000.0" => 41.22
                    )
                ),
                "formula" => "H2O"
            ),
            "CO2" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 37.13,
                        "500.0" => 44.63,
                        "1000.0" => 54.31
                    )
                ),
                "formula" => "CO2"
            ),
            "O3" => Dict(
                "properties" => Dict(
                    "heat_capacity" => Dict(
                        "298.15" => 39.20,
                        "500.0" => 45.04,
                        "1000.0" => 53.35
                    )
                ),
                "formula" => "O3"
            )
        )
        
        # Save as JSON
        json_str = JSON.json(sample_data, 4)
        prepare_sample_data("tde", "tde_data.json", json_str)
        println("TDE sample data created.")
        return true
    end
end

"""
    prepare_thermoml_data()

Prepare ThermoML sample data.
"""
function prepare_thermoml_data()
    println("Preparing ThermoML data...")
    
    # Copy sample file
    if copy_sample_file("thermoml", "sample_thermoml.dat", "thermoml.xml")
        println("ThermoML data prepared successfully.")
        return true
    else
        # Create a minimal sample XML file
        sample_data = """
        <?xml version="1.0" encoding="UTF-8"?>
        <DataReport xmlns="http://www.iupac.org/namespaces/ThermoML">
          <Version>2.0</Version>
          <Citation>
            <Title>ThermoML Sample Data</Title>
            <Year>2023</Year>
          </Citation>
          <Compound>
            <ID>1</ID>
            <CASRegistry>7727-37-9</CASRegistry>
            <RegNum>N2</RegNum>
            <Name>Nitrogen</Name>
            <Formula>N2</Formula>
          </Compound>
          <Compound>
            <ID>2</ID>
            <CASRegistry>7782-44-7</CASRegistry>
            <RegNum>O2</RegNum>
            <Name>Oxygen</Name>
            <Formula>O2</Formula>
          </Compound>
          <Compound>
            <ID>3</ID>
            <CASRegistry>7732-18-5</CASRegistry>
            <RegNum>H2O</RegNum>
            <Name>Water</Name>
            <Formula>H2O</Formula>
          </Compound>
          <PureOrMixtureData>
            <ID>1</ID>
            <PropertyGroup>
              <Property>
                <nPropNumber>1</nPropNumber>
                <PropName>Heat capacity at constant pressure</PropName>
                <Measurement>
                  <ID>1</ID>
                  <Method>
                    <Name>DSC</Name>
                  </Method>
                  <Temperature>
                    <Value>298.15</Value>
                  </Temperature>
                  <Value>29.12</Value>
                  <Uncertainty>0.1</Uncertainty>
                </Measurement>
              </Property>
            </PropertyGroup>
          </PureOrMixtureData>
        </DataReport>
        """
        
        prepare_sample_data("thermoml", "thermoml.xml", sample_data)
        println("ThermoML sample data created.")
        return true
    end
end

"""
    prepare_janaf_data()

Prepare JANAF Thermochemical Tables sample data.
"""
function prepare_janaf_data()
    println("Preparing JANAF Thermochemical Tables data...")
    
    # Copy sample file
    if copy_sample_file("janaf", "sample_janaf.dat", "data.txt")
        println("JANAF data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = """
        JANAF THERMOCHEMICAL TABLES
        NITROGEN (N2)
        T/K    Cp/J mol-1 K-1    S/J mol-1 K-1    H/kJ mol-1
        298.15    29.124    191.61    0.0
        300    29.165    191.79    0.054
        400    29.246    198.99    3.167
        500    29.580    204.64    6.086
        600    30.110    209.52    9.063
        700    30.754    213.74    12.122
        800    31.433    217.45    15.237
        900    32.068    220.72    18.444
        1000    32.698    223.76    21.749
        
        OXYGEN (O2)
        T/K    Cp/J mol-1 K-1    S/J mol-1 K-1    H/kJ mol-1
        298.15    29.376    205.15    0.0
        300    29.390    205.33    0.054
        400    30.105    213.87    3.067
        500    31.091    220.97    6.085
        600    32.090    226.90    9.185
        700    32.981    232.04    12.369
        800    33.733    236.53    15.629
        900    34.409    240.62    18.848
        1000    34.870    244.30    22.542
        
        WATER (H2O)
        T/K    Cp/J mol-1 K-1    S/J mol-1 K-1    H/kJ mol-1
        298.15    33.590    188.84    -241.826
        300    33.596    189.01    -241.782
        400    34.262    198.79    -238.403
        500    35.226    207.42    -238.918
        600    36.325    215.01    -232.149
        700    37.495    221.76    -229.201
        800    38.721    227.85    -226.077
        900    39.987    233.45    -222.776
        1000    41.217    238.65    -219.300
        """
        
        prepare_sample_data("janaf", "data.txt", sample_data)
        println("JANAF sample data created.")
        return true
    end
end

"""
    prepare_nasa_cea_data()

Prepare NASA CEA Database sample data.
"""
function prepare_nasa_cea_data()
    println("Preparing NASA CEA Database data...")
    
    # Copy sample file
    if copy_sample_file("nasa-cea", "sample_nasa7.dat", "thermo.inp")
        println("NASA CEA data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = """
        thermo
         N2               N 2                 G  200.000  6000.000  1000.00     1
          2.9526600E+00  1.3968500E-03 -4.9262400E-07  7.8600100E-11 -4.6074700E-15    2
         -9.2395500E+02  5.8718000E+00  3.5310000E+00 -1.2366600E-04 -5.0299900E-07    3
          2.4353000E-09 -1.4088100E-12 -1.0469800E+03  2.9674700E+00                   4
         O2               O 2                 G  200.000  6000.000  1000.00     1
          3.6609600E+00  6.5636500E-04 -1.4114900E-07  2.0579800E-11 -1.2991300E-15    2
         -1.2159800E+03  3.4153600E+00  3.7824600E+00 -2.9967300E-03  9.8473000E-06    3
         -9.6813000E-09  3.2437300E-12 -1.0639400E+03  3.6576800E+00                   4
         H2O               H 2  O 1            G  200.000  6000.000  1000.00     1
          2.6770400E+00  2.9731800E-03 -7.7376900E-07  9.4433700E-11 -4.2690100E-15    2
         -2.9886000E+04  6.8825600E+00  4.1986400E+00 -2.0364000E-03  6.5203400E-06    3
         -5.4879300E-09  1.7719700E-12 -3.0293700E+04 -8.4900900E-01                   4
         CO2               C 1  O 2            G  200.000  6000.000  1000.00     1
          4.6365100E+00  2.7414600E-03 -9.9589800E-07  1.6038700E-10 -9.1619900E-15    2
         -4.9024900E+04 -1.9349000E+00  2.3567800E+00  8.9846000E-03 -7.1235600E-06    3
          2.4591900E-09 -1.4370000E-13 -4.8372000E+04  9.9010500E+00                   4
         O3                O 3                 G  200.000  6000.000  1000.00     1
          1.2330300E+01 -1.1932500E-02  7.9874100E-06 -1.7719500E-09  1.2607600E-13    2
          1.2675600E+04 -4.0882300E+01  3.4073800E+00  2.0537900E-03 -1.3828300E-06    3
         -6.2799400E-10  6.5435400E-13  1.5864500E+04  8.2815200E+00                   4
        END
        """
        
        prepare_sample_data("nasa-cea", "thermo.inp", sample_data)
        println("NASA CEA sample data created.")
        return true
    end
end

"""
    prepare_chemkin_data()

Prepare CHEMKIN format sample data.
"""
function prepare_chemkin_data()
    println("Preparing CHEMKIN format data...")
    
    # Copy sample file
    if copy_sample_file("chemkin", "sample_chemkin.dat", "therm.dat")
        println("CHEMKIN data prepared successfully.")
        return true
    else
        # Create a minimal sample
        sample_data = """
        THERMO
           300.000  1000.000  5000.000
        ! CHEMKIN format thermodynamic data
        N2                      N   2               G   200.000  6000.000 1000.00      1
         2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
        -9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
         2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00                   4
        O2                      O   2               G   200.000  6000.000 1000.00      1
         3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
        -1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
        -9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00                   4
        H2O                     H   2O   1          G   200.000  6000.000 1000.00      1
         2.67703787E+00 2.97318329E-03-7.73769690E-07 9.44336689E-11-4.26900959E-15    2
        -2.98858938E+04 6.88255571E+00 4.19863520E+00-2.03640170E-03 6.52034160E-06    3
        -5.48792690E-09 1.77196800E-12-3.02937260E+04-8.49009010E-01                   4
        CO2                     C   1O   2          G   200.000  6000.000 1000.00      1
         4.63651110E+00 2.74145690E-03-9.95897590E-07 1.60386660E-10-9.16198570E-15    2
        -4.90249040E+04-1.93489550E+00 2.35677520E+00 8.98459680E-03-7.12356270E-06    3
         2.45919220E-09-1.43699550E-13-4.83719710E+04 9.90105220E+00                   4
        O3                      O   3               G   200.000  6000.000 1000.00      1
         1.23302914E+01-1.19324783E-02 7.98741278E-06-1.77194552E-09 1.26075824E-13    2
         1.26755831E+04-4.08823374E+01 3.40738221E+00 2.05379063E-03-1.38283480E-06    3
        -6.27993680E-10 6.54353945E-13 1.58644979E+04 8.28152238E+00                   4
        END
        """
        
        prepare_sample_data("chemkin", "therm.dat", sample_data)
        println("CHEMKIN sample data created.")
        return true
    end
end

"""
    prepare_gri_mech_data()

Prepare GRI-MECH 3.0 sample data.
"""
function prepare_gri_mech_data()
    println("Preparing GRI-MECH 3.0 data...")
    
    # Copy sample file if exists in cache first
    existing_file = joinpath(dirname(dirname(@__FILE__)), "data", "cache", "gri-mech", "therm.dat")
    if isfile(existing_file)
        println("  Using existing GRI-MECH data file.")
        return true
    end
    
    # Check external directory
    external_file = joinpath(dirname(dirname(@__FILE__)), "data", "external", "chemkin", "therm.dat")
    if isfile(external_file)
        # Create target directory
        cache_dir = joinpath(dirname(dirname(@__FILE__)), "data", "cache", "gri-mech")
        mkpath(cache_dir)
        
        # Copy file
        target_file = joinpath(cache_dir, "therm.dat")
        cp(external_file, target_file, force=true)
        
        println("  Copied GRI-MECH data from external directory.")
        return true
    end
    
    # Create a minimal sample (similar to CHEMKIN format)
    sample_data = """
    THERMO
       300.000  1000.000  5000.000
    ! GRI-MECH 3.0 thermodynamic data
    N2                      N   2               G   200.000  6000.000 1000.00      1
     2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
    -9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
     2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00                   4
    O2                      O   2               G   200.000  6000.000 1000.00      1
     3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
    -1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
    -9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00                   4
    H2O                     H   2O   1          G   200.000  6000.000 1000.00      1
     2.67703787E+00 2.97318329E-03-7.73769690E-07 9.44336689E-11-4.26900959E-15    2
    -2.98858938E+04 6.88255571E+00 4.19863520E+00-2.03640170E-03 6.52034160E-06    3
    -5.48792690E-09 1.77196800E-12-3.02937260E+04-8.49009010E-01                   4
    CO2                     C   1O   2          G   200.000  6000.000 1000.00      1
     4.63651110E+00 2.74145690E-03-9.95897590E-07 1.60386660E-10-9.16198570E-15    2
    -4.90249040E+04-1.93489550E+00 2.35677520E+00 8.98459680E-03-7.12356270E-06    3
     2.45919220E-09-1.43699550E-13-4.83719710E+04 9.90105220E+00                   4
    CH4                     C   1H   4          G   200.000  6000.000 1000.00      1
     1.63552643E+00 1.00842795E-02-3.36916254E-06 5.34958667E-10-3.15518833E-14    2
    -1.00125611E+04 9.99371680E+00 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
    -4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
    END
    """
    
    prepare_sample_data("gri-mech", "therm.dat", sample_data)
    println("GRI-MECH sample data created.")
    return true
end

"""
    main()

Main function to prepare test data for all sources.
"""
function main()
    println("JThermodynamicsData - Test Data Preparation")
    println("===========================================")
    println("Started at: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    
    # Create cache directories
    create_cache_directories()
    
    # Prepare sample data for each source
    sources = [
        ("ATcT", prepare_atct_data),
        ("Burcat", prepare_burcat_data),
        ("NIST WebBook", prepare_nist_webbook_data),
        ("TDE", prepare_tde_data),
        ("ThermoML", prepare_thermoml_data),
        ("JANAF", prepare_janaf_data),
        ("NASA CEA", prepare_nasa_cea_data),
        ("CHEMKIN", prepare_chemkin_data),
        ("GRI-MECH", prepare_gri_mech_data)
    ]
    
    success_count = 0
    
    for (name, prepare_fn) in sources
        println("\n$(repeat("=", 50))")
        if prepare_fn()
            success_count += 1
        end
    end
    
    println("\nData preparation summary:")
    println("  Successfully prepared $(success_count) out of $(length(sources)) data sources.")
    
    println("\nNow that all test data is prepared, you can run the full workflow:")
    println("  1. julia scripts/fetch_all_sources.jl")
    println("  2. julia run_all_species_plots.jl")
    println("  3. julia scripts/sync_json_to_database.jl")
    println("\nOr run everything in one command:")
    println("  julia run_full_workflow.jl")
end

# Run the main function
main()