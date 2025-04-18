"""
Parsers for NASA 7-coefficient and NASA 9-coefficient polynomials.
"""

"""
    parse_nasa7_polynomial(lines::Vector{String})

Parse NASA 7-coefficient polynomial format.
NASA 7-coefficient format uses two temperature ranges with 7 coefficients each.
"""
function parse_nasa7_polynomial(lines::Vector{String})
    if length(lines) < 3
        error("Invalid NASA-7 format: not enough lines")
    end
    
    # Extract species information from first line
    header = lines[1]
    species_name = strip(header[1:18])
    
    # Parse temperature ranges from second line
    temp_line = lines[2]
    temp_ranges = []
    
    if length(temp_line) >= 45
        try
            t_low = parse(Float64, strip(temp_line[1:10]))
            t_mid = parse(Float64, strip(temp_line[11:20]))
            t_high = parse(Float64, strip(temp_line[21:30]))
            
            push!(temp_ranges, (t_low, t_mid))
            push!(temp_ranges, (t_mid, t_high))
        catch e
            error("Failed to parse temperature ranges: $(e)")
        end
    else
        error("Invalid temperature range line format")
    end
    
    # Parse coefficients
    coefficients = []
    
    # High-temperature coefficients (first set in NASA format)
    high_coefs = zeros(7)
    if length(lines) >= 4
        # Read first line of coefficients
        if length(lines[3]) >= 75
            high_coefs[1] = parse(Float64, strip(lines[3][1:15]))
            high_coefs[2] = parse(Float64, strip(lines[3][16:30]))
            high_coefs[3] = parse(Float64, strip(lines[3][31:45]))
            high_coefs[4] = parse(Float64, strip(lines[3][46:60]))
            high_coefs[5] = parse(Float64, strip(lines[3][61:75]))
        end
        
        # Read second line of coefficients
        if length(lines[4]) >= 45
            high_coefs[6] = parse(Float64, strip(lines[4][1:15]))
            high_coefs[7] = parse(Float64, strip(lines[4][16:30]))
        end
    end
    
    # Low-temperature coefficients (second set in NASA format)
    low_coefs = zeros(7)
    if length(lines) >= 6
        # Read first line of coefficients
        if length(lines[5]) >= 75
            low_coefs[1] = parse(Float64, strip(lines[5][1:15]))
            low_coefs[2] = parse(Float64, strip(lines[5][16:30]))
            low_coefs[3] = parse(Float64, strip(lines[5][31:45]))
            low_coefs[4] = parse(Float64, strip(lines[5][46:60]))
            low_coefs[5] = parse(Float64, strip(lines[5][61:75]))
        end
        
        # Read second line of coefficients
        if length(lines[6]) >= 45
            low_coefs[6] = parse(Float64, strip(lines[6][1:15]))
            low_coefs[7] = parse(Float64, strip(lines[6][16:30]))
        end
    end
    
    # In NASA format, the high-temp coeffs are listed first, but they're for the second range
    push!(coefficients, low_coefs)
    push!(coefficients, high_coefs)
    
    # Extract molecular formula if available
    formula = ""
    if length(header) >= 45
        formula = strip(header[25:44])
    end
    
    # Extract other metadata
    metadata = Dict{String, Any}()
    if length(header) >= 75
        metadata["phase"] = strip(header[45:45])
        metadata["source"] = strip(header[57:64])
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

"""
    parse_nasa9_polynomial(lines::Vector{String})

Parse NASA 9-coefficient polynomial format.
NASA 9-coefficient format uses multiple temperature ranges with 9 coefficients each.
"""
function parse_nasa9_polynomial(lines::Vector{String})
    if length(lines) < 4  # Minimum lines for a single range
        error("Invalid NASA-9 format: not enough lines")
    end
    
    # Extract species information from first line
    header = lines[1]
    species_name = ""
    formula = ""
    
    if length(header) >= 24
        species_name = strip(header[1:24])
    end
    
    if length(header) >= 80
        formula = strip(header[25:80])
    end
    
    # Second line contains number of temperature intervals and additional info
    if length(lines[2]) < 2
        error("Invalid NASA-9 format: second line too short")
    end
    
    n_intervals = 0
    try
        n_intervals = parse(Int, strip(lines[2][1:2]))
    catch e
        error("Failed to parse number of temperature intervals: $(e)")
    end
    
    if n_intervals <= 0
        error("Invalid number of temperature intervals: $n_intervals")
    end
    
    # Parse common temperature points (n_intervals + 1 points)
    temps = []
    temp_line = lines[3]
    
    for i in 0:n_intervals
        start_idx = 1 + i * 11
        end_idx = start_idx + 10
        
        if length(temp_line) >= end_idx
            try
                temp = parse(Float64, strip(temp_line[start_idx:end_idx]))
                push!(temps, temp)
            catch e
                error("Failed to parse temperature point $i: $(e)")
            end
        else
            error("Temperature line not long enough for all intervals")
        end
    end
    
    # Create temperature ranges
    temp_ranges = [(temps[i], temps[i+1]) for i in 1:length(temps)-1]
    
    # Parse coefficients for each interval
    coefficients = []
    line_idx = 4
    
    for interval in 1:n_intervals
        if line_idx + 2 > length(lines)
            error("Not enough lines for interval $interval coefficients")
        end
        
        # NASA-9 uses 3 lines per interval with 3 coefficients per line
        coefs = zeros(9)
        
        for coef_line in 0:2
            for coef_idx in 0:2
                if coef_line * 3 + coef_idx >= 9
                    break
                end
                
                start_idx = 1 + coef_idx * 16
                end_idx = start_idx + 15
                
                if length(lines[line_idx + coef_line]) >= end_idx
                    coefs[coef_line * 3 + coef_idx + 1] = parse(Float64, 
                                                          strip(lines[line_idx + coef_line][start_idx:end_idx]))
                end
            end
        end
        
        push!(coefficients, coefs)
        line_idx += 3
    end
    
    return Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature_ranges" => temp_ranges,
        "coefficients" => coefficients,
        "polynomial_type" => "nasa9",
        "metadata" => Dict{String, Any}()
    )
end

"""
    calculate_nasa7_cp(coeffs::Vector{Float64}, T::Float64)

Calculate specific heat capacity (Cp/R) using NASA 7-coefficient polynomial.
"""
function calculate_nasa7_cp(coeffs::Vector{Float64}, T::Float64)
    return coeffs[1] + coeffs[2]*T + coeffs[3]*T^2 + coeffs[4]*T^3 + coeffs[5]*T^4
end

"""
    calculate_nasa7_enthalpy(coeffs::Vector{Float64}, T::Float64)

Calculate enthalpy (H/RT) using NASA 7-coefficient polynomial.
"""
function calculate_nasa7_enthalpy(coeffs::Vector{Float64}, T::Float64)
    return coeffs[1] + coeffs[2]*T/2 + coeffs[3]*T^2/3 + coeffs[4]*T^3/4 + coeffs[5]*T^4/5 + coeffs[6]/T
end

"""
    calculate_nasa7_entropy(coeffs::Vector{Float64}, T::Float64)

Calculate entropy (S/R) using NASA 7-coefficient polynomial.
"""
function calculate_nasa7_entropy(coeffs::Vector{Float64}, T::Float64)
    return coeffs[1]*log(T) + coeffs[2]*T + coeffs[3]*T^2/2 + coeffs[4]*T^3/3 + coeffs[5]*T^4/4 + coeffs[7]
end

"""
    calculate_nasa9_cp(coeffs::Vector{Float64}, T::Float64)

Calculate specific heat capacity (Cp/R) using NASA 9-coefficient polynomial.
"""
function calculate_nasa9_cp(coeffs::Vector{Float64}, T::Float64)
    return coeffs[1]*T^(-2) + coeffs[2]*T^(-1) + coeffs[3] + coeffs[4]*T + coeffs[5]*T^2 + 
           coeffs[6]*T^3 + coeffs[7]*T^4
end

"""
    calculate_nasa9_enthalpy(coeffs::Vector{Float64}, T::Float64)

Calculate enthalpy (H/RT) using NASA 9-coefficient polynomial.
"""
function calculate_nasa9_enthalpy(coeffs::Vector{Float64}, T::Float64)
    return -coeffs[1]*T^(-2) + coeffs[2]*log(T)/T + coeffs[3] + coeffs[4]*T/2 + 
           coeffs[5]*T^2/3 + coeffs[6]*T^3/4 + coeffs[7]*T^4/5 + coeffs[8]/T
end

"""
    calculate_nasa9_entropy(coeffs::Vector{Float64}, T::Float64)

Calculate entropy (S/R) using NASA 9-coefficient polynomial.
"""
function calculate_nasa9_entropy(coeffs::Vector{Float64}, T::Float64)
    return -coeffs[1]*T^(-2)/2 - coeffs[2]*T^(-1) + coeffs[3]*log(T) + coeffs[4]*T + 
           coeffs[5]*T^2/2 + coeffs[6]*T^3/3 + coeffs[7]*T^4/4 + coeffs[9]
end

"""
    nasa_polynomial_to_thermodynamic_data(data::Dict, source::String, reliability_score::Float64=0.0)

Convert parsed NASA polynomial data to ThermodynamicData structure.
"""
function nasa_polynomial_to_thermodynamic_data(data::Dict, source::String, reliability_score::Float64=0.0)
    species_name = data["species_name"]
    formula = data["formula"]
    temp_ranges = data["temperature_ranges"]
    coefficients = data["coefficients"]
    poly_type = data["polynomial_type"]
    
    # Create PolynomialCoefficients structures
    polynomial_type = poly_type == "nasa7" ? NASA7 : NASA9
    
    polynomials = []
    for i in 1:length(temp_ranges)
        poly = PolynomialCoefficients(
            polynomial_type,
            [temp_ranges[i]],
            [coefficients[i]],
            load_data_source(source),
            reliability_score,
            string(Dates.now())
        )
        push!(polynomials, poly)
    end
    
    return ThermodynamicData(
        species_name,
        formula,
        get(data["metadata"], "cas", ""),
        get(data["metadata"], "molecular_weight", 0.0),
        polynomials,
        nothing,
        0,  # Source priority to be set later
        NONE,  # Uncertainty method
        Dict{String, Any}(),
        data["metadata"]
    )
end

"""
    extract_formula_elements(formula::String)

Extract chemical elements and their counts from a formula string.
"""
function extract_formula_elements(formula::String)
    elements = Dict{String, Int}()
    
    # Simple regex for formula parsing
    # This is a basic implementation that might need enhancement for complex formulas
    pattern = r"([A-Z][a-z]*)(\d*)"
    matches = eachmatch(pattern, formula)
    
    for m in matches
        element = m.captures[1]
        count = m.captures[2]
        
        if isempty(count)
            count = "1"
        end
        
        elements[element] = get(elements, element, 0) + parse(Int, count)
    end
    
    return elements
end

"""
    calculate_molecular_weight(formula::String)

Calculate molecular weight from a chemical formula.
"""
function calculate_molecular_weight(formula::String)
    elements = extract_formula_elements(formula)
    
    mw = 0.0
    for (element, count) in elements
        if haskey(ATOMIC_MASSES, element)
            mw += ATOMIC_MASSES[element] * count
        else
            @warn "Unknown element in formula: $element"
        end
    end
    
    return mw
end

"""
    parse_nasa7_file(file_path::String)

Parse a NASA 7-coefficient thermodynamic data file and return a dictionary of species data.
"""
function parse_nasa7_file(file_path::String)
    if !isfile(file_path)
        error("NASA file not found: $file_path")
    end
    
    lines = readlines(file_path)
    
    # Dictionary to store species data
    species_data = Dict{String, Dict}()
    
    # NASA 7-coefficient files typically have 4 lines per species
    # Example:
    # N2               N 2    G100.0    5000.0  1000.0      1
    #  2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
    #  9.89224864E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
    #  2.43530612E-09-1.40881235E-12 9.46834564E+02 2.96747038E+00                   4
    
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        
        # Skip empty lines or comments
        if isempty(line) || startswith(line, "!")
            i += 1
            continue
        end
        
        # Try to parse a record (4 lines per species)
        if i + 3 <= length(lines)
            species_lines = lines[i:i+3]
            
            # Try to parse the species data
            try
                species_data_dict = parse_nasa7_polynomial(species_lines)
                
                # Add to the species data dictionary
                species_name = species_data_dict["species_name"]
                species_data[species_name] = species_data_dict
                
                # Move to the next species
                i += 4
            catch e
                # If parsing fails, try moving to the next line
                @warn "Failed to parse NASA 7 polynomial at line $i: $e"
                i += 1
            end
        else
            # Not enough lines remaining for a full record
            break
        end
    end
    
    return species_data
end

"""
    nasa_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a NASA thermodynamic data file and store the data in the database.
"""
function nasa_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_nasa7_file(file_path)
    
    count = 0
    
    # Process data within a transaction for efficiency
    transaction(function()
        for (species_name, data) in species_data
            # Get molecular weight from metadata or calculate from formula
            mw = get(data["metadata"], "molecular_weight", 0.0)
            if mw == 0.0 && !isempty(data["formula"])
                mw = calculate_molecular_weight(data["formula"])
            end
            
            # Get CAS number from metadata
            cas = get(data["metadata"], "cas", "")
            
            # Add or get species in the database
            species_id = add_species(conn, species_name, data["formula"], cas, mw, data["metadata"])
            
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
            
            count += 1
        end
    end, conn)
    
    return count
end

"""
    parse_nasa9_file(file_path::String)

Parse a NASA 9-coefficient thermodynamic data file and return a dictionary of species data.
"""
function parse_nasa9_file(file_path::String)
    if !isfile(file_path)
        error("NASA file not found: $file_path")
    end
    
    lines = readlines(file_path)
    
    # Dictionary to store species data
    species_data = Dict{String, Dict}()
    
    # NASA 9-coefficient files are more complex and have variable number of lines per species
    # The first line has the species name and formula
    # The second line has the number of temperature intervals
    # The third line has the temperature points
    # Then, 3 lines per interval for the 9 coefficients (3 per line)
    
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        
        # Skip empty lines or comments
        if isempty(line) || startswith(line, "!")
            i += 1
            continue
        end
        
        # Try to parse the species data
        try
            # Check if we have at least 4 lines (minimum for a NASA 9 record with 1 interval)
            if i + 3 <= length(lines)
                # First line: species name and formula
                header = lines[i]
                
                # Second line: number of temperature intervals
                intervals_line = lines[i+1]
                n_intervals = 0
                
                try
                    n_intervals = parse(Int, strip(intervals_line[1:2]))
                catch e
                    @warn "Failed to parse number of temperature intervals at line $(i+1)"
                    i += 1
                    continue
                end
                
                if n_intervals <= 0
                    @warn "Invalid number of temperature intervals: $n_intervals at line $(i+1)"
                    i += 1
                    continue
                end
                
                # Check if we have enough lines for all intervals
                if i + 2 + 3*n_intervals <= length(lines)
                    # Extract all lines for this species
                    species_lines = lines[i:i+2+3*n_intervals]
                    
                    # Parse the NASA 9 polynomial
                    species_data_dict = parse_nasa9_polynomial(species_lines)
                    
                    # Add to the species data dictionary
                    species_name = species_data_dict["species_name"]
                    species_data[species_name] = species_data_dict
                    
                    # Move past this species record
                    i += 3 + 3*n_intervals
                else
                    @warn "Not enough lines for all intervals at line $i"
                    i += 1
                end
            else
                i += 1
            end
        catch e
            # If parsing fails, try moving to the next line
            @warn "Failed to parse NASA 9 polynomial at line $i: $e"
            i += 1
        end
    end
    
    return species_data
end

"""
    nasa9_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a NASA 9-coefficient thermodynamic data file and store the data in the database.
"""
function nasa9_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_nasa9_file(file_path)
    
    count = 0
    
    # Process data within a transaction for efficiency
    transaction(function()
        for (species_name, data) in species_data
            # Get molecular weight from metadata or calculate from formula
            mw = get(data["metadata"], "molecular_weight", 0.0)
            if mw == 0.0 && !isempty(data["formula"])
                mw = calculate_molecular_weight(data["formula"])
            end
            
            # Get CAS number from metadata
            cas = get(data["metadata"], "cas", "")
            
            # Add or get species in the database
            species_id = add_species(conn, species_name, data["formula"], cas, mw, data["metadata"])
            
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
                "nasa9", 
                temp_min, 
                temp_max, 
                thermo_data,
                Dict(),  # No uncertainty data
                reliability_score
            )
            
            count += 1
        end
    end, conn)
    
    return count
end