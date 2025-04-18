"""
Parser for JANAF Thermochemical Tables.
JANAF tables provide tabular thermodynamic data at standard temperature points.
"""

"""
    parse_janaf_file(file_path::String)

Parse a JANAF thermochemical table file and return a dictionary of species data.
"""
function parse_janaf_file(file_path::String)
    if !isfile(file_path)
        error("JANAF file not found: $file_path")
    end
    
    lines = readlines(file_path)
    return parse_janaf_data(lines)
end

"""
    parse_janaf_data(lines::Vector{String})

Parse JANAF thermochemical data from a vector of lines.
"""
function parse_janaf_data(lines::Vector{String})
    # JANAF tables have a specific format with header information
    # and tabular data at standard temperature points
    if length(lines) < 10
        error("Invalid JANAF file: not enough lines")
    end
    
    # Extract species information from header
    species_name = ""
    formula = ""
    phase = ""
    cas_number = ""
    
    # Parse header lines to extract metadata
    for i in 1:10
        if i > length(lines)
            break
        end
        
        line = lines[i]
        
        if occursin("Formula:", line) || occursin("FORMULA:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                formula = strip(parts[2])
            end
        elseif occursin("Name:", line) || occursin("NAME:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                species_name = strip(parts[2])
            end
        elseif occursin("State:", line) || occursin("STATE:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                phase = strip(parts[2])
            end
        elseif occursin("CAS:", line) || occursin("CAS Number:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                cas_number = strip(parts[2])
            end
        end
    end
    
    # Find the start of the data table
    table_start = 0
    for i in 1:length(lines)
        if occursin("T(K)", lines[i]) || occursin("Temperature", lines[i])
            table_start = i
            break
        end
    end
    
    if table_start == 0
        error("Could not find data table in JANAF file")
    end
    
    # Parse column headers
    header_line = lines[table_start]
    columns = []
    
    if occursin("T(K)", header_line)
        push!(columns, "T")
    end
    if occursin("Cp", header_line)
        push!(columns, "Cp")
    end
    if occursin("S", header_line) && !occursin("ΔfS", header_line)
        push!(columns, "S")
    end
    if occursin("-(G-H(Tr))/T", header_line) || occursin("-(G⁰-H⁰(Tr))/T", header_line)
        push!(columns, "-(G-H)/T")
    end
    if occursin("H-H(Tr)", header_line) || occursin("H⁰-H⁰(Tr)", header_line)
        push!(columns, "H-H(Tr)")
    end
    if occursin("ΔfH", header_line)
        push!(columns, "ΔfH")
    end
    if occursin("ΔfG", header_line)
        push!(columns, "ΔfG")
    end
    if occursin("log Kf", header_line)
        push!(columns, "log_Kf")
    end
    
    # Parse data table
    data_table = Dict{String, Vector{Float64}}()
    for col in columns
        data_table[col] = Float64[]
    end
    
    for i in (table_start + 1):length(lines)
        line = lines[i]
        
        # Skip empty lines or lines with non-numeric data
        if isempty(strip(line)) || !any(isdigit, line)
            continue
        end
        
        # Parse data line
        # JANAF tables typically use fixed-width columns
        # This is a simplified parser that assumes consistent column widths
        values = split(line)
        
        if length(values) >= length(columns)
            for j in 1:length(columns)
                col = columns[j]
                try
                    val = parse(Float64, values[j])
                    push!(data_table[col], val)
                catch e
                    # Skip non-numeric values
                    @warn "Failed to parse value in column $col at line $i: $(values[j])"
                end
            end
        end
    end
    
    # Check if we have temperature data
    if !haskey(data_table, "T") || isempty(data_table["T"])
        error("No temperature data found in JANAF table")
    end
    
    # Metadata
    metadata = Dict{String, Any}(
        "cas" => cas_number,
        "phase" => phase,
        "source" => "JANAF"
    )
    
    # Get temperature range
    temps = data_table["T"]
    temp_min = minimum(temps)
    temp_max = maximum(temps)
    
    # Calculate molecular weight
    mw = 0.0
    if !isempty(formula)
        mw = calculate_molecular_weight(formula)
        metadata["molecular_weight"] = mw
    end
    
    return Dict(
        "species_name" => species_name,
        "formula" => formula,
        "tabular_data" => data_table,
        "temperature_range" => [temp_min, temp_max],
        "metadata" => metadata
    )
end

"""
    janaf_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a JANAF table file and store the data in the database.
"""
function janaf_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    data = parse_janaf_file(file_path)
    
    # Use transaction for database operations
    transaction(function()
        # Get or create species
        species_name = data["species_name"]
        formula = data["formula"]
        cas = get(data["metadata"], "cas", "")
        mw = get(data["metadata"], "molecular_weight", 0.0)
        
        species_id = add_species(conn, species_name, formula, cas, mw, data["metadata"])
    
    # Add thermodynamic data
    temp_range = data["temperature_range"]
    temp_min = temp_range[1]
    temp_max = temp_range[2]
    
    # Create data dictionary
    thermo_data = Dict(
        "tabular_data" => data["tabular_data"]
    )
    
    # Add to database
        add_thermodynamic_data(
            conn, 
            species_id, 
            source_name, 
            "tabular", 
            temp_min, 
            temp_max, 
            thermo_data,
            Dict(),  # No uncertainty data
            reliability_score
        )
        return species_name
    end, conn)
    
    return data["species_name"]
end

"""
    janaf_to_nasa7(data::Dict)

Convert JANAF tabular data to NASA 7-coefficient polynomials.
"""
function janaf_to_nasa7(data::Dict)
    tabular_data = data["tabular_data"]
    
    # Check if we have the required data
    if !haskey(tabular_data, "T") || !haskey(tabular_data, "Cp") || 
       !haskey(tabular_data, "H-H(Tr)") || !haskey(tabular_data, "S")
        error("Missing required data columns in JANAF table")
    end
    
    # Extract data
    temps = tabular_data["T"]
    cps = tabular_data["Cp"]
    hs = tabular_data["H-H(Tr)"]
    ss = tabular_data["S"]
    
    # Sort data by temperature
    indices = sortperm(temps)
    temps = temps[indices]
    cps = cps[indices]
    hs = hs[indices]
    ss = ss[indices]
    
    # Split into two temperature ranges
    # Typically, 200-1000K and 1000-6000K for NASA polynomials
    split_temp = 1000.0
    
    # Find index closest to 1000K
    split_idx = findmin(abs.(temps .- split_temp))[2]
    
    # Low temperature range
    low_temps = temps[1:split_idx]
    low_cps = cps[1:split_idx]
    low_hs = hs[1:split_idx]
    low_ss = ss[1:split_idx]
    
    # High temperature range
    high_temps = temps[split_idx:end]
    high_cps = cps[split_idx:end]
    high_hs = hs[split_idx:end]
    high_ss = ss[split_idx:end]
    
    # Fit NASA 7-coefficient polynomials using least squares
    low_coeffs = fit_nasa7_coefficients(low_temps, low_cps, low_hs, low_ss)
    high_coeffs = fit_nasa7_coefficients(high_temps, high_cps, high_hs, high_ss)
    
    # Create temperature ranges
    temp_ranges = [(minimum(low_temps), maximum(low_temps)), (minimum(high_temps), maximum(high_temps))]
    
    # Return NASA 7 polynomial data
    return Dict(
        "species_name" => data["species_name"],
        "formula" => data["formula"],
        "temperature_ranges" => temp_ranges,
        "coefficients" => [low_coeffs, high_coeffs],
        "polynomial_type" => "nasa7",
        "metadata" => data["metadata"]
    )
end

"""
    fit_nasa7_coefficients(temps::Vector{Float64}, cps::Vector{Float64}, 
                         hs::Vector{Float64}, ss::Vector{Float64})

Fit NASA 7-coefficient polynomials to JANAF tabular data.
"""
function fit_nasa7_coefficients(temps::Vector{Float64}, cps::Vector{Float64}, 
                              hs::Vector{Float64}, ss::Vector{Float64})
    # NASA 7-coefficient polynomials:
    # Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    # H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    # S/R = a1*log(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    
    # This is a simplified fitting approach
    # In practice, more sophisticated fitting techniques are used
    
    # Fit Cp/R polynomial (a1-a5)
    n = length(temps)
    A = zeros(n, 5)
    
    for i in 1:n
        t = temps[i]
        A[i, 1] = 1.0
        A[i, 2] = t
        A[i, 3] = t^2
        A[i, 4] = t^3
        A[i, 5] = t^4
    end
    
    # Convert Cp to Cp/R
    cp_r = cps ./ R
    
    # Solve least squares problem for a1-a5
    a1_a5 = A \ cp_r
    
    # Fit a6 using H/RT
    h_rt = hs ./ (R .* temps)
    h_computed = zeros(n)
    
    for i in 1:n
        t = temps[i]
        h_computed[i] = a1_a5[1] + a1_a5[2]*t/2 + a1_a5[3]*t^2/3 + a1_a5[4]*t^3/4 + a1_a5[5]*t^4/5
    end
    
    # Solve for a6
    a6 = sum((h_rt - h_computed) .* temps) / n
    
    # Fit a7 using S/R
    s_r = ss ./ R
    s_computed = zeros(n)
    
    for i in 1:n
        t = temps[i]
        s_computed[i] = a1_a5[1]*log(t) + a1_a5[2]*t + a1_a5[3]*t^2/2 + a1_a5[4]*t^3/3 + a1_a5[5]*t^4/4
    end
    
    # Solve for a7
    a7 = sum(s_r - s_computed) / n
    
    # Return coefficients
    return [a1_a5..., a6, a7]
end