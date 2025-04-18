"""
Parser for the NIST ThermoData Engine (TDE) format data.
TDE provides evaluated thermodynamic property data with uncertainty information.
"""

"""
    parse_tde_file(file_path::String)

Parse a TDE export file and return a dictionary of thermodynamic data.
"""
function parse_tde_file(file_path::String)
    if !isfile(file_path)
        error("TDE file not found: $file_path")
    end
    
    # TDE files can be in different formats depending on the export options
    # We'll attempt to handle common formats like XML, JSON, and ThermoML
    
    # Check file extension
    if endswith(lowercase(file_path), ".xml")
        # Might be ThermoML format
        try
            return parse_thermoml_file(file_path)
        catch e
            @warn "Failed to parse as ThermoML, trying generic XML: $e"
            
            # Try to parse as generic XML
            return parse_tde_xml(file_path)
        end
    elseif endswith(lowercase(file_path), ".json")
        # JSON format
        return parse_tde_json(file_path)
    else
        # Try to determine format from content
        content = read(file_path, String)
        
        if startswith(strip(content), "<")
            # Likely XML format
            try
                return parse_thermoml_file(file_path)
            catch e
                @warn "Failed to parse as ThermoML, trying generic XML: $e"
                
                # Try to parse as generic XML
                return parse_tde_xml(file_path)
            end
        elseif startswith(strip(content), "{") || startswith(strip(content), "[")
            # Likely JSON format
            return parse_tde_json(file_path)
        else
            # Try to parse as a tabular format
            return parse_tde_tabular(file_path)
        end
    end
end

"""
    parse_tde_xml(file_path::String)

Parse a TDE XML file format.
"""
function parse_tde_xml(file_path::String)
    # Parse XML document
    xdoc = LightXML.parse_file(file_path)
    root = LightXML.root(xdoc)
    
    species_data = Dict{String, Dict}()
    
    # Extract compound information
    for compound_elem in LightXML.get_elements_by_tagname(root, "Compound")
        compound_data = Dict{String, Any}()
        
        # Extract compound identifier
        species_name = ""
        formula = ""
        cas_number = ""
        
        # Try different tag names for compound name
        for name_tag in ["Name", "CanonicalName", "CompoundName", "StandardName", "PreferredName"]
            name_elem = LightXML.find_element(compound_elem, name_tag)
            if name_elem !== nothing
                species_name = LightXML.content(name_elem)
                if !isempty(species_name)
                    break
                end
            end
        end
        
        # Try different tag names for formula
        for formula_tag in ["Formula", "MolecularFormula", "ChemicalFormula"]
            formula_elem = LightXML.find_element(compound_elem, formula_tag)
            if formula_elem !== nothing
                formula = LightXML.content(formula_elem)
                if !isempty(formula)
                    break
                end
            end
        end
        
        # Try different tag names for CAS number
        for cas_tag in ["CASNo", "CASNumber", "CASRegistryNumber", "CAS"]
            cas_elem = LightXML.find_element(compound_elem, cas_tag)
            if cas_elem !== nothing
                cas_number = LightXML.content(cas_elem)
                if !isempty(cas_number)
                    break
                end
            end
        end
        
        compound_data["species_name"] = species_name
        compound_data["formula"] = formula
        compound_data["metadata"] = Dict{String, Any}(
            "cas" => cas_number,
            "source" => "TDE"
        )
        
        # Extract property data
        tabular_data = Dict{String, Vector{Float64}}(
            "T" => Float64[]
        )
        
        # Find all property elements for this compound
        for prop_elem in LightXML.get_elements_by_tagname(compound_elem, "Property")
            # Get property name
            prop_name = ""
            name_elem = LightXML.find_element(prop_elem, "Name")
            if name_elem !== nothing
                prop_name = LightXML.content(name_elem)
            else
                # Try attribute
                if LightXML.has_attribute(prop_elem, "name")
                    prop_name = LightXML.attribute(prop_elem, "name")
                end
            end
            
            # Skip if no property name
            if isempty(prop_name)
                continue
            end
            
            # Map property name to standard key
            prop_key = ""
            if occursin("heat capacity", lowercase(prop_name)) || occursin("cp", lowercase(prop_name))
                prop_key = "Cp"
            elseif occursin("enthalpy", lowercase(prop_name)) || occursin("h", lowercase(prop_name))
                prop_key = "H"
            elseif occursin("entropy", lowercase(prop_name)) || occursin("s", lowercase(prop_name))
                prop_key = "S"
            elseif occursin("gibbs", lowercase(prop_name)) || occursin("g", lowercase(prop_name))
                prop_key = "G"
            else
                # Skip properties we're not interested in
                continue
            end
            
            # Initialize vectors for this property
            if !haskey(tabular_data, prop_key)
                tabular_data[prop_key] = Float64[]
                tabular_data["$(prop_key)_uncertainty"] = Float64[]
            end
            
            # Find data points
            for point_elem in LightXML.get_elements_by_tagname(prop_elem, "DataPoint")
                # Get temperature
                temp = 0.0
                temp_elem = LightXML.find_element(point_elem, "Temperature")
                if temp_elem !== nothing
                    val_elem = LightXML.find_element(temp_elem, "Value")
                    if val_elem !== nothing
                        try
                            temp = parse(Float64, LightXML.content(val_elem))
                        catch e
                            @warn "Failed to parse temperature: $e"
                            continue
                        end
                    end
                end
                
                # Skip if no temperature data
                if temp == 0.0
                    continue
                end
                
                # Get property value
                value = 0.0
                value_elem = LightXML.find_element(point_elem, "Value")
                if value_elem !== nothing
                    try
                        value = parse(Float64, LightXML.content(value_elem))
                    catch e
                        @warn "Failed to parse property value: $e"
                        continue
                    end
                end
                
                # Get uncertainty if available
                uncertainty = 0.0
                uncert_elem = LightXML.find_element(point_elem, "Uncertainty")
                if uncert_elem !== nothing
                    try
                        uncertainty = parse(Float64, LightXML.content(uncert_elem))
                    catch e
                        @warn "Failed to parse uncertainty: $e"
                    end
                end
                
                # Add temperature if not already in the list
                if !(temp in tabular_data["T"])
                    push!(tabular_data["T"], temp)
                    
                    # Add placeholder values for other properties at this temperature
                    for p in keys(tabular_data)
                        if p != "T" && p != prop_key && p != "$(prop_key)_uncertainty"
                            push!(tabular_data[p], NaN)
                        end
                    end
                end
                
                # Find index of this temperature
                idx = findfirst(x -> x == temp, tabular_data["T"])
                
                # Add property value
                if length(tabular_data[prop_key]) < idx
                    # Pad with NaNs if needed
                    append!(tabular_data[prop_key], fill(NaN, idx - length(tabular_data[prop_key])))
                end
                
                tabular_data[prop_key][idx] = value
                
                # Add uncertainty if available
                if uncertainty > 0.0
                    if length(tabular_data["$(prop_key)_uncertainty"]) < idx
                        append!(tabular_data["$(prop_key)_uncertainty"], 
                               fill(NaN, idx - length(tabular_data["$(prop_key)_uncertainty"])))
                    end
                    
                    tabular_data["$(prop_key)_uncertainty"][idx] = uncertainty
                end
            end
        end
        
        # Add tabular data if we found any
        if !isempty(tabular_data["T"])
            # Get temperature range
            temps = tabular_data["T"]
            temp_min = minimum(temps)
            temp_max = maximum(temps)
            
            compound_data["tabular_data"] = tabular_data
            compound_data["temperature_range"] = [temp_min, temp_max]
            
            # Calculate molecular weight if formula is available
            if !isempty(formula)
                mw = calculate_molecular_weight(formula)
                compound_data["metadata"]["molecular_weight"] = mw
            end
            
            # Add to species data dictionary
            species_data[species_name] = compound_data
        end
    end
    
    # Free XML document
    LightXML.free(xdoc)
    
    return species_data
end

"""
    parse_tde_json(file_path::String)

Parse a TDE JSON file format.
"""
function parse_tde_json(file_path::String)
    # Parse JSON file
    json_data = JSON.parsefile(file_path)
    
    species_data = Dict{String, Dict}()
    
    # Determine JSON structure format
    if isa(json_data, Dict) && haskey(json_data, "compounds")
        # Process compounds array
        compounds = json_data["compounds"]
        
        for compound in compounds
            species_name = get(compound, "name", "")
            formula = get(compound, "formula", "")
            cas_number = get(compound, "cas", "")
            
            if isempty(species_name)
                continue
            end
            
            compound_data = Dict{String, Any}(
                "species_name" => species_name,
                "formula" => formula,
                "metadata" => Dict{String, Any}(
                    "cas" => cas_number,
                    "source" => "TDE"
                )
            )
            
            # Process properties if available
            if haskey(compound, "properties")
                properties = compound["properties"]
                
                tabular_data = Dict{String, Vector{Float64}}(
                    "T" => Float64[]
                )
                
                # Process property data
                for (prop_name, prop_data) in properties
                    # Map property name to standard key
                    prop_key = ""
                    if occursin("heat capacity", lowercase(prop_name)) || occursin("cp", lowercase(prop_name))
                        prop_key = "Cp"
                    elseif occursin("enthalpy", lowercase(prop_name)) || occursin("h", lowercase(prop_name))
                        prop_key = "H"
                    elseif occursin("entropy", lowercase(prop_name)) || occursin("s", lowercase(prop_name))
                        prop_key = "S"
                    elseif occursin("gibbs", lowercase(prop_name)) || occursin("g", lowercase(prop_name))
                        prop_key = "G"
                    else
                        # Skip properties we're not interested in
                        continue
                    end
                    
                    # Initialize vectors for this property
                    if !haskey(tabular_data, prop_key)
                        tabular_data[prop_key] = Float64[]
                        tabular_data["$(prop_key)_uncertainty"] = Float64[]
                    end
                    
                    # Process data points
                    if isa(prop_data, Vector)
                        for point in prop_data
                            # Get temperature
                            temp = get(point, "temperature", 0.0)
                            
                            # Skip if no temperature data
                            if temp == 0.0
                                continue
                            end
                            
                            # Get property value
                            value = get(point, "value", 0.0)
                            
                            # Get uncertainty if available
                            uncertainty = get(point, "uncertainty", 0.0)
                            
                            # Add temperature if not already in the list
                            if !(temp in tabular_data["T"])
                                push!(tabular_data["T"], temp)
                                
                                # Add placeholder values for other properties at this temperature
                                for p in keys(tabular_data)
                                    if p != "T" && p != prop_key && p != "$(prop_key)_uncertainty"
                                        push!(tabular_data[p], NaN)
                                    end
                                end
                            end
                            
                            # Find index of this temperature
                            idx = findfirst(x -> x == temp, tabular_data["T"])
                            
                            # Add property value
                            if length(tabular_data[prop_key]) < idx
                                # Pad with NaNs if needed
                                append!(tabular_data[prop_key], fill(NaN, idx - length(tabular_data[prop_key])))
                            end
                            
                            tabular_data[prop_key][idx] = value
                            
                            # Add uncertainty if available
                            if uncertainty > 0.0
                                if length(tabular_data["$(prop_key)_uncertainty"]) < idx
                                    append!(tabular_data["$(prop_key)_uncertainty"], 
                                           fill(NaN, idx - length(tabular_data["$(prop_key)_uncertainty"])))
                                end
                                
                                tabular_data["$(prop_key)_uncertainty"][idx] = uncertainty
                            end
                        end
                    end
                end
                
                # Add tabular data if we found any
                if !isempty(tabular_data["T"])
                    # Get temperature range
                    temps = tabular_data["T"]
                    temp_min = minimum(temps)
                    temp_max = maximum(temps)
                    
                    compound_data["tabular_data"] = tabular_data
                    compound_data["temperature_range"] = [temp_min, temp_max]
                    
                    # Calculate molecular weight if formula is available
                    if !isempty(formula)
                        mw = calculate_molecular_weight(formula)
                        compound_data["metadata"]["molecular_weight"] = mw
                    end
                    
                    # Add to species data dictionary
                    species_data[species_name] = compound_data
                end
            end
        end
    end
    
    return species_data
end

"""
    parse_tde_tabular(file_path::String)

Parse a TDE tabular text file format.
"""
function parse_tde_tabular(file_path::String)
    lines = readlines(file_path)
    
    species_data = Dict{String, Dict}()
    
    # Try to identify the file structure
    # TDE tabular exports often have a header section and data section
    
    # Extract species information from header lines
    species_name = ""
    formula = ""
    cas_number = ""
    
    for i in 1:min(20, length(lines))
        line = lines[i]
        
        if occursin("Compound:", line) || occursin("Name:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                species_name = strip(parts[2])
            end
        elseif occursin("Formula:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                formula = strip(parts[2])
            end
        elseif occursin("CAS RN:", line) || occursin("CAS:", line)
            parts = split(line, ":")
            if length(parts) >= 2
                cas_number = strip(parts[2])
            end
        end
    end
    
    # Find the start of the data table
    table_start = 0
    for i in 1:length(lines)
        if occursin("T/K", lines[i]) || occursin("Temperature", lines[i])
            table_start = i
            break
        end
    end
    
    if table_start == 0 || isempty(species_name)
        # Could not identify file structure
        @warn "Could not identify TDE tabular file structure"
        return species_data
    end
    
    # Parse column headers
    header_line = lines[table_start]
    headers = split(header_line)
    
    # Map headers to standard property keys
    header_map = Dict{String, String}()
    
    for header in headers
        if occursin("T/K", header) || occursin("Temperature", header)
            header_map[header] = "T"
        elseif occursin("Cp", header) || occursin("Heat Capacity", header)
            header_map[header] = "Cp"
        elseif occursin("H", header) && !occursin("Error", header)
            header_map[header] = "H"
        elseif occursin("S", header) && !occursin("Error", header)
            header_map[header] = "S"
        elseif occursin("G", header) && !occursin("Error", header)
            header_map[header] = "G"
        elseif occursin("U(Cp)", header) || occursin("Error(Cp)", header)
            header_map[header] = "Cp_uncertainty"
        elseif occursin("U(H)", header) || occursin("Error(H)", header)
            header_map[header] = "H_uncertainty"
        elseif occursin("U(S)", header) || occursin("Error(S)", header)
            header_map[header] = "S_uncertainty"
        elseif occursin("U(G)", header) || occursin("Error(G)", header)
            header_map[header] = "G_uncertainty"
        end
    end
    
    # Check if we have temperature column
    if !any(values(header_map) .== "T")
        @warn "No temperature column found in TDE tabular file"
        return species_data
    end
    
    # Initialize tabular data
    tabular_data = Dict{String, Vector{Float64}}()
    for mapped_header in unique(values(header_map))
        tabular_data[mapped_header] = Float64[]
    end
    
    # Parse data lines
    for i in (table_start + 1):length(lines)
        line = strip(lines[i])
        
        if isempty(line) || !any(isdigit, line)
            continue
        end
        
        values = split(line)
        
        if length(values) >= length(headers)
            for j in 1:length(headers)
                header = headers[j]
                
                if haskey(header_map, header)
                    mapped_header = header_map[header]
                    
                    try
                        value = parse(Float64, values[j])
                        push!(tabular_data[mapped_header], value)
                    catch e
                        @warn "Failed to parse value in column $header at line $i: $(values[j])"
                        push!(tabular_data[mapped_header], NaN)
                    end
                end
            end
        end
    end
    
    # Check if we have any data
    if isempty(tabular_data["T"])
        @warn "No data found in TDE tabular file"
        return species_data
    end
    
    # Create compound data
    compound_data = Dict{String, Any}(
        "species_name" => species_name,
        "formula" => formula,
        "metadata" => Dict{String, Any}(
            "cas" => cas_number,
            "source" => "TDE"
        )
    )
    
    # Get temperature range
    temps = tabular_data["T"]
    temp_min = minimum(temps)
    temp_max = maximum(temps)
    
    compound_data["tabular_data"] = tabular_data
    compound_data["temperature_range"] = [temp_min, temp_max]
    
    # Calculate molecular weight if formula is available
    if !isempty(formula)
        mw = calculate_molecular_weight(formula)
        compound_data["metadata"]["molecular_weight"] = mw
    end
    
    # Add to species data dictionary
    species_data[species_name] = compound_data
    
    return species_data
end

"""
    tde_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a TDE data file and store the data in the database.
"""
function tde_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_tde_file(file_path)
    
    count = 0
    
    # Transaction for bulk inserts
    transaction(function()
        for (species_name, data) in species_data
            # Skip if no tabular data
            if !haskey(data, "tabular_data") || isempty(data["tabular_data"]) || !haskey(data["tabular_data"], "T")
                continue
            end
            
            # Get or create species
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
            
            # Create uncertainty dictionary
            uncertainty = Dict{String, Any}()
            
            # Extract uncertainty data
            for key in keys(data["tabular_data"])
                if endswith(key, "_uncertainty")
                    base_key = key[1:end-12]
                    uncertainty[base_key] = data["tabular_data"][key]
                end
            end
            
            # Add to database
            add_thermodynamic_data(
                conn, 
                species_id, 
                source_name, 
                "tabular", 
                temp_min, 
                temp_max, 
                thermo_data,
                uncertainty,
                reliability_score
            )
            
            count += 1
        end
    end, conn)
    
    return count
end