"""
    JThermodynamicsData.parsers.nist

This module handles parsing thermodynamic data from the NIST WebBook.
It retrieves and processes data from the NIST Chemistry WebBook to extract
thermodynamic properties and convert them to NASA-7 polynomial format.
"""

module NISTParser

using HTTP
using Gumbo
using Cascadia
using Dates
using Printf
using Statistics

# Export the main function
export parse_nist_webbook

"""
    parse_nist_webbook(cas_number::String; debug::Bool=false)

Parse thermodynamic data for a chemical species from NIST WebBook given its CAS number.
Returns a dictionary with NASA-7 polynomial coefficients and metadata.

Args:
    cas_number: CAS Registry Number for the species
    debug: Enable debug output
"""
function parse_nist_webbook(cas_number::String; debug::Bool=false)
    # Format the URL for NIST WebBook
    url = "https://webbook.nist.gov/cgi/cbook.cgi?ID=$(cas_number)&Units=SI&Mask=1"
    
    if debug
        println("Fetching NIST WebBook data from: $url")
    end
    
    # Fetch the HTML content
    response = HTTP.get(url)
    html = String(response.body)
    parsed_html = parsehtml(html)
    
    # Extract the chemical name - fixed to avoid Cascadia.findfirst issue
    species_name = extract_chemical_name_fixed(html)
    if debug
        println("Extracted species name: $species_name")
    end
    
    # Extract thermodynamic data from tables
    tables = extract_tables_fixed(html)
    
    # Initialize data containers
    temperature_points = Float64[]
    cp_values = Float64[]
    h_values = Float64[]
    s_values = Float64[]
    
    # Find the heat capacity table (Shomate Equation coefficients table)
    shomate_table = find_shomate_table(tables)
    
    if shomate_table === nothing
        if debug
            println("No Shomate equation table found")
        end
        throw(ErrorException("Shomate equation table not found in NIST WebBook for CAS $cas_number"))
    end
    
    if debug
        println("Found Shomate equation table")
    end
    
    # Parse the Shomate coefficients
    shomate_coeffs = parse_shomate_coefficients(shomate_table)
    
    if isempty(shomate_coeffs)
        if debug
            println("Failed to extract Shomate coefficients")
        end
        throw(ErrorException("Failed to extract Shomate coefficients from NIST WebBook for CAS $cas_number"))
    end
    
    if debug
        println("Extracted Shomate coefficients: $(length(shomate_coeffs)) temperature ranges")
        for (temp_range, coeffs) in shomate_coeffs
            println("Temperature range: $(temp_range[1]) - $(temp_range[2]) K")
            println("Coefficients: $coeffs")
        end
    end
    
    # Extract standard enthalpy of formation
    h_formation = extract_formation_enthalpy_fixed(html, debug)
    if debug
        println("Extracted formation enthalpy: $h_formation kJ/mol")
    end
    
    # Generate temperature points and calculate thermodynamic properties
    if !isempty(shomate_coeffs)
        # Get the full temperature range
        min_temp = minimum([range[1] for (range, _) in shomate_coeffs])
        max_temp = maximum([range[2] for (range, _) in shomate_coeffs])
        
        # Generate temperature points
        step = (max_temp - min_temp) / 50
        temperature_points = collect(min_temp:step:max_temp)
        
        # Calculate Cp, H, and S at each temperature point
        for T in temperature_points
            # Find the appropriate temperature range and coefficients
            range_idx = findfirst(x -> x[1][1] <= T <= x[1][2], shomate_coeffs)
            if range_idx !== nothing
                _, coeffs = shomate_coeffs[range_idx]
                
                # Calculate values using Shomate equation
                t = T / 1000.0  # Shomate equations use t = T/1000
                
                # Heat capacity: Cp = A + B*t + C*t^2 + D*t^3 + E/t^2
                cp = coeffs[1] + coeffs[2]*t + coeffs[3]*t^2 + coeffs[4]*t^3 + coeffs[5]/(t^2)
                
                # Enthalpy: H-H298 = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F - H
                h = coeffs[1]*t + coeffs[2]*t^2/2 + coeffs[3]*t^3/3 + coeffs[4]*t^4/4 - 
                    coeffs[5]/t + coeffs[6] - coeffs[8]
                
                # Entropy: S = A*ln(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G
                s = coeffs[1]*log(t) + coeffs[2]*t + coeffs[3]*t^2/2 + coeffs[4]*t^3/3 - 
                    coeffs[5]/(2*t^2) + coeffs[7]
                
                push!(cp_values, cp)
                push!(h_values, h)
                push!(s_values, s)
            end
        end
    end
    
    # If we couldn't extract or calculate data points, throw an error
    if isempty(temperature_points) || isempty(cp_values)
        if debug
            println("Failed to generate thermodynamic data points")
        end
        throw(ErrorException("Failed to generate thermodynamic data for CAS $cas_number"))
    end
    
    if debug
        println("Generated $(length(temperature_points)) temperature points")
        println("Temperature range: $(minimum(temperature_points)) - $(maximum(temperature_points)) K")
    end
    
    # Create tabular data format
    tabular_data = Dict(
        "temperatures" => temperature_points,
        "heat_capacities" => cp_values,
        "enthalpies" => h_values, 
        "entropies" => s_values,
        "reference_enthalpy" => h_formation
    )
    
    # Define the temperature ranges for NASA polynomials
    # Use the original Shomate ranges if available
    if !isempty(shomate_coeffs)
        if length(shomate_coeffs) >= 2
            temp_ranges = [shomate_coeffs[1][1], shomate_coeffs[2][1]]
        else
            # Handle single temperature range case
            temp_range = shomate_coeffs[1][1]
            midpoint = (temp_range[1] + temp_range[2]) / 2
            temp_ranges = [[temp_range[1], midpoint], [midpoint, temp_range[2]]]
        end
    else
        # Otherwise use standard ranges or adapt to the data
        temp_min = minimum(temperature_points)
        temp_max = maximum(temperature_points)
        
        if temp_min <= 300.0 && temp_max >= 5000.0
            # Standard temperature ranges
            temp_ranges = [[200.0, 1000.0], [1000.0, 6000.0]]
        else
            # Adapt to available data
            midpoint = (temp_min + temp_max) / 2
            temp_ranges = [[temp_min, midpoint], [midpoint, temp_max]]
        end
    end
    
    # Return structured data
    result = Dict(
        "source" => "NIST-WEBBOOK",
        "priority" => 11,
        "reliability_score" => 9.0,
        "polynomial_type" => "nasa7",
        "temperature_ranges" => temp_ranges,
        "tabular_data" => tabular_data,
        "shomate_coefficients" => shomate_coeffs,
        "species_name" => species_name,
        "formula" => species_name,  # Use name as placeholder
        "uncertainty" => 0.03  # 3% - typical for NIST data
    )
    
    if debug
        println("Extracted NIST WebBook data for $species_name")
    end
    
    return result
end

"""
    extract_chemical_name_fixed(html_str)

Extract the chemical name from the NIST WebBook HTML string.
Uses regex instead of DOM parsing to avoid issues with Cascadia.
"""
function extract_chemical_name_fixed(html_str)
    # Try to find the title tag content
    title_match = match(r"<title>\s*(.*?)\s*</title>", html_str)
    if title_match !== nothing
        title_text = title_match.captures[1]
        # Format: "Species Name - NIST Chemistry WebBook"
        if occursin(" - NIST", title_text)
            return strip(split(title_text, " - ")[1])
        end
        return strip(title_text)
    end
    
    # Fallback: look for the first h1 element
    h1_match = match(r"<h1[^>]*>\s*(.*?)\s*</h1>", html_str)
    if h1_match !== nothing
        return strip(h1_match.captures[1])
    end
    
    # Last resort: return "Unknown Species"
    return "Unknown Species"
end

"""
    extract_tables_fixed(html_str)

Extract all tables from the HTML string.
Returns a list of table HTML strings.
"""
function extract_tables_fixed(html_str)
    tables = []
    # Find all table tags and their content
    table_pattern = r"<table[^>]*>(.*?)</table>"s  # 's' flag for dot to match newlines
    for m in eachmatch(table_pattern, html_str)
        push!(tables, m.match)
    end
    return tables
end

"""
    find_shomate_table(tables)

Find the Shomate equation coefficients table from the list of tables.
"""
function find_shomate_table(tables)
    for table in tables
        # Check if this is the Shomate table by looking for specific headers or content
        if occursin("Shomate", table) && (
            occursin("Temperature (K)", table) || 
            occursin("A", table) && occursin("B", table) && occursin("C", table)
           )
            return table
        end
    end
    return nothing
end

"""
    parse_shomate_coefficients(table_html)

Parse the Shomate equation coefficients from the table HTML.
Returns a list of (temperature_range, coefficients) tuples.
"""
function parse_shomate_coefficients(table_html)
    result = []
    
    # Extract temperature ranges
    temp_ranges = []
    temp_pattern = r">\s*(\d+)\.?\s+to\s+(\d+)\.?\s*<"
    for m in eachmatch(temp_pattern, table_html)
        if length(m.captures) >= 2
            low = parse(Float64, m.captures[1])
            high = parse(Float64, m.captures[2])
            push!(temp_ranges, [low, high])
        end
    end
    
    # Extract coefficients A through H
    coeff_patterns = [
        r">\s*A\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*B\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*C\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*D\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*E\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*F\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*G\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>",
        r">\s*H\s*</th>\s*<td[^>]*>\s*(-?\d+\.\d+)\s*</td>"
    ]
    
    # Process each temperature range
    for i in 1:length(temp_ranges)
        coeffs = []
        
        # Starting position for searching this column's data
        start_pos = 1
        for pattern in coeff_patterns
            m = match(pattern, table_html, start_pos)
            if m !== nothing
                push!(coeffs, parse(Float64, m.captures[1]))
                start_pos = m.offset + length(m.match)  # Move past this match
            else
                # If we can't find a coefficient, use a default value
                push!(coeffs, 0.0)
            end
        end
        
        # Only add if we found some coefficients
        if !isempty(coeffs)
            push!(result, (temp_ranges[i], coeffs))
        end
    end
    
    return result
end

"""
    extract_formation_enthalpy_fixed(html_str, debug=false)

Extract the standard enthalpy of formation from the NIST WebBook HTML string.
Uses regex instead of DOM parsing to avoid issues with Cascadia.
"""
function extract_formation_enthalpy_fixed(html_str, debug=false)
    # Look for tables with enthalpy of formation
    enthalpy_patterns = [
        r"<td[^>]*>\s*Δ\s*<sub>\s*f\s*</sub>\s*H°\s*<sub>\s*gas\s*</sub>.*?<td[^>]*>\s*(-?[\d\.]+)\s*</td>",
        r"<td[^>]*>\s*Standard enthalpy of formation.*?<td[^>]*>\s*(-?[\d\.]+)\s*</td>",
        r"<td[^>]*>\s*Heat of formation.*?<td[^>]*>\s*(-?[\d\.]+)\s*</td>"
    ]
    
    for pattern in enthalpy_patterns
        m = match(pattern, html_str)
        if m !== nothing
            value = parse(Float64, m.captures[1])
            
            # Check units (kJ/mol vs kcal/mol)
            units_check = match(r"(kJ/mol|kcal/mol)", html_str, m.offset)
            if units_check !== nothing && units_check.match == "kcal/mol"
                value *= 4.184  # Convert kcal to kJ
            end
            
            if debug
                println("Found enthalpy of formation: $value kJ/mol")
            end
            
            return value
        end
    end
    
    # If still not found, look for it in text paragraphs
    paragraph_pattern = r"enthalpy of formation[^<>]*?(-?[\d\.]+)\s*(kJ/mol|kcal/mol)"
    m = match(paragraph_pattern, html_str)
    if m !== nothing
        value = parse(Float64, m.captures[1])
        
        # Convert if needed
        if m.captures[2] == "kcal/mol"
            value *= 4.184  # Convert kcal to kJ
        end
        
        if debug
            println("Found enthalpy of formation in paragraph: $value kJ/mol")
        end
        
        return value
    end
    
    # If still not found, use a sensible default
    if debug
        println("Enthalpy of formation not found, using default value")
    end
    return 0.0  # Default value when not found
end

end # module