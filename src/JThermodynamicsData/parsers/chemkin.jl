"""
Parser for CHEMKIN-format thermodynamic data.
"""

"""
    parse_chemkin_file(file_path::String)

Parse a CHEMKIN-format thermodynamic data file and return a dictionary of species data.
"""
function parse_chemkin_file(file_path::String)
    if !isfile(file_path)
        error("CHEMKIN file not found: $file_path")
    end
    
    lines = readlines(file_path)
    return parse_chemkin_data(lines)
end

"""
    parse_chemkin_data(lines::Vector{String})

Parse CHEMKIN-format thermodynamic data from a vector of lines.
"""
function parse_chemkin_data(lines::Vector{String})
    species_data = Dict{String, Dict}()
    
    # Skip header lines until THERMO or THERMO ALL
    i = 1
    header_end = 0
    while i <= length(lines)
        if occursin(r"^THERMO\s*(ALL)?", uppercase(strip(lines[i])))
            header_end = i
            break
        end
        i += 1
    end
    
    if header_end == 0
        error("THERMO section not found in CHEMKIN file")
    end
    
    # Check for temperature ranges in the line after THERMO
    default_temp_ranges = []
    if header_end + 1 <= length(lines)
        temp_line = lines[header_end + 1]
        if length(temp_line) >= 30 && all(isdigit, strip(temp_line[1:5]))
            try
                t_low = parse(Float64, strip(temp_line[1:10]))
                t_mid = parse(Float64, strip(temp_line[11:20]))
                t_high = parse(Float64, strip(temp_line[21:30]))
                
                push!(default_temp_ranges, (t_low, t_mid))
                push!(default_temp_ranges, (t_mid, t_high))
                
                # Skip the temperature line
                i = header_end + 2
            catch e
                # Not a valid temperature line, continue parsing
                i = header_end + 1
            end
        else
            # No temperature line, continue parsing
            i = header_end + 1
        end
    else
        # No more lines after THERMO
        return species_data
    end
    
    # Parse species data
    while i <= length(lines)
        # Check for END marker
        if occursin(r"^END", uppercase(strip(lines[i])))
            break
        end
        
        # Need at least 4 lines for a species (header + 3 coefficient lines)
        if i + 3 > length(lines)
            break
        end
        
        # Read 4 lines per species
        species_lines = lines[i:min(i+3, length(lines))]
        
        if length(species_lines) < 4
            i += length(species_lines)
            continue
        end
        
        # Parse the species using NASA 7-coefficient format
        try
            data = parse_chemkin_nasa7_polynomial(species_lines)
            species_name = data["species_name"]
            
            # If no temperature ranges were in the data, use the default
            if isempty(data["temperature_ranges"]) && !isempty(default_temp_ranges)
                data["temperature_ranges"] = default_temp_ranges
            end
            
            # Add to species data dictionary
            species_data[species_name] = data
        catch e
            @warn "Failed to parse species at line $i: $e"
        end
        
        i += 4
    end
    
    return species_data
end

"""
    chemkin_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a CHEMKIN file and store the data in the database.
"""
function chemkin_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_chemkin_file(file_path)
    
    # Use transaction for database operations
    transaction(function()
        for (species_name, data) in species_data
            # Get or create species
            molecular_weight = calculate_molecular_weight(data["formula"])
            species_id = add_species(conn, species_name, data["formula"], "", molecular_weight)
            
            # Add thermodynamic data
            temp_ranges = data["temperature_ranges"]
            
            if isempty(temp_ranges)
                @warn "No temperature ranges found for species: $species_name"
                continue
            end
            
            # Temperature range for database (min to max across all ranges)
            temp_min = minimum([range[1] for range in temp_ranges])
            temp_max = maximum([range[2] for range in temp_ranges])
            
            # Create data dictionary
            thermo_data = Dict(
                "coefficients" => data["coefficients"],
                "temperature_ranges" => temp_ranges
            )
            
            # Add to database
            add_thermodynamic_data(
                conn, 
                species_id, 
                source_name, 
                "nasa7", 
                temp_min, 
                temp_max, 
                thermo_data,
                Dict(),  # No uncertainty data
                reliability_score
            )
        end
    end, conn)
    
    return length(species_data)
end

"""
    parse_chemkin_nasa7_polynomial(species_lines::Vector{String})

Parse a species entry in CHEMKIN-format NASA 7-coefficient polynomial.
"""
function parse_chemkin_nasa7_polynomial(species_lines::Vector{String})
    if length(species_lines) < 4
        error("Invalid CHEMKIN entry: not enough lines")
    end
    
    # First line contains species name and metadata
    header = species_lines[1]
    species_name = strip(header[1:16])  # CHEMKIN format puts name in first 16 chars
    
    # Parse temperature ranges
    t_low = 0.0
    t_mid = 0.0
    t_high = 0.0
    
    # Get temperature ranges from header line
    if length(header) >= 80
        try
            t_min = parse(Float64, strip(header[46:55]))
            t_max = parse(Float64, strip(header[56:65]))
            t_mid = parse(Float64, strip(header[66:75]))
            
            # CHEMKIN specific format handling - sometimes the order is different
            if t_mid > t_max
                # If the "midpoint" is higher than max, then order is min, max, mid
                t_low = t_min
                t_high = t_max
                t_mid = t_mid
            else
                # Otherwise standard order min, mid, max
                t_low = t_min
                t_high = t_max
                t_mid = t_mid
            end
        catch e
            # If we can't parse temperatures from the header, use defaults
            t_low = 200.0
            t_mid = 1000.0
            t_high = 5000.0
        end
    else
        # Use defaults if header is too short
        t_low = 200.0
        t_mid = 1000.0
        t_high = 5000.0
    end
    
    # Get molecular formula
    formula = ""
    if length(header) >= 45
        # Extract formula from elements specification
        elements = []
        formula_start = 24
        
        # Scan fields in 5-character blocks (CHEMKIN element + count format)
        for i in 0:3
            pos = formula_start + i*5
            if pos + 4 <= length(header)
                element = strip(header[pos:pos+1])
                count_str = strip(header[pos+2:pos+4])
                
                if !isempty(element) && element != "00"
                    count = 0
                    try
                        count = parse(Int, count_str)
                    catch
                        count = 1
                    end
                    
                    if count > 0
                        push!(elements, (element, count))
                    end
                end
            end
        end
        
        # Build formula string
        formula = join([e * (c > 1 ? string(c) : "") for (e, c) in elements], "")
    end
    
    # Parse coefficients - CHEMKIN format has high-temp coeffs on lines 2-3, low-temp on lines 3-4
    # High-temperature coefficients (first two lines)
    high_coefs = zeros(7)
    
    # Extract 5 coefficients from line 2
    line2 = species_lines[2]
    if length(line2) >= 75
        for i in 0:4
            pos = 1 + i*15
            high_coefs[i+1] = parse(Float64, strip(line2[pos:pos+14]))
        end
    end
    
    # Extract 2 coefficients from line 3 (first 2 fields)
    line3 = species_lines[3]
    if length(line3) >= 30
        for i in 0:1
            pos = 1 + i*15
            high_coefs[i+6] = parse(Float64, strip(line3[pos:pos+14]))
        end
    end
    
    # Low-temperature coefficients (last line and part of previous)
    low_coefs = zeros(7)
    
    # Extract 5 coefficients from line 3 (last 3 fields) and line 4 (first 2 fields)
    if length(line3) >= 75
        for i in 0:2
            pos = 31 + i*15
            if pos + 14 <= length(line3)
                low_coefs[i+1] = parse(Float64, strip(line3[pos:pos+14]))
            end
        end
    end
    
    # Get last 4 coefficients from line 4
    line4 = species_lines[4]
    if length(line4) >= 60
        for i in 0:3
            pos = 1 + i*15
            if i < 2
                # First 2 coefficients (index 4-5)
                low_coefs[i+4] = parse(Float64, strip(line4[pos:pos+14]))
            else
                # Last 2 coefficients (index 6-7)
                low_coefs[i+2] = parse(Float64, strip(line4[pos:pos+14]))
            end
        end
    end
    
    # Create temperature ranges
    temp_ranges = [(t_low, t_mid), (t_mid, t_high)]
    
    # Package coefficients - high temp first in CHEMKIN format
    coefficients = [high_coefs, low_coefs]
    
    # Extract phase and other metadata if available
    metadata = Dict{String, Any}()
    if length(header) >= 45
        phase = header[45:45]
        if !isempty(phase) && phase != " "
            metadata["phase"] = phase
        end
    end
    
    return Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature_ranges" => temp_ranges,
        "coefficients" => coefficients,
        "polynomial_type" => "nasa7",
        "metadata" => metadata
    )
end