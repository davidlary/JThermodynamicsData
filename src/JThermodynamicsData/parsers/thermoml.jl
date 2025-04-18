"""
Parser for ThermoML format data.
ThermoML is an XML-based IUPAC standard format for thermodynamic data.
"""

"""
    parse_thermoml_file(file_path::String)

Parse a ThermoML XML file and return a dictionary of thermodynamic data.
"""
function parse_thermoml_file(file_path::String)
    if !isfile(file_path)
        error("ThermoML file not found: $file_path")
    end
    
    # Parse XML document
    xdoc = LightXML.parse_file(file_path)
    root = LightXML.root(xdoc)
    
    # Check if this is a valid ThermoML document
    if LightXML.name(root) != "ThermoML"
        error("Not a valid ThermoML document")
    end
    
    # Extract data
    data = parse_thermoml_document(root)
    
    # Free XML document
    LightXML.free(xdoc)
    
    return data
end

"""
    parse_thermoml_document(root::LightXML.XMLElement)

Parse a ThermoML XML document and extract thermodynamic data.
"""
function parse_thermoml_document(root::LightXML.XMLElement)
    species_data = Dict{String, Dict}()
    
    # Get version information
    version = ""
    if LightXML.has_attribute(root, "version")
        version = LightXML.attribute(root, "version")
    end
    
    # Extract citation information
    citation = Dict{String, Any}()
    citation_elem = LightXML.find_element(root, "Citation")
    
    if citation_elem !== nothing
        title_elem = LightXML.find_element(citation_elem, "Title")
        if title_elem !== nothing
            citation["title"] = LightXML.content(title_elem)
        end
        
        authors = []
        for author_elem in LightXML.get_elements_by_tagname(citation_elem, "Author")
            name = ""
            
            first_name = LightXML.find_element(author_elem, "FirstName")
            if first_name !== nothing
                name = LightXML.content(first_name)
            end
            
            last_name = LightXML.find_element(author_elem, "LastName")
            if last_name !== nothing
                if !isempty(name)
                    name *= " "
                end
                name *= LightXML.content(last_name)
            end
            
            if !isempty(name)
                push!(authors, name)
            end
        end
        
        citation["authors"] = authors
        
        year_elem = LightXML.find_element(citation_elem, "Year")
        if year_elem !== nothing
            citation["year"] = LightXML.content(year_elem)
        end
        
        source_elem = LightXML.find_element(citation_elem, "Source")
        if source_elem !== nothing
            citation["source"] = LightXML.content(source_elem)
        end
    end
    
    # Extract compound information
    for compound_elem in LightXML.get_elements_by_tagname(root, "Compound")
        compound_data = Dict{String, Any}()
        
        # Extract compound identifier
        reg_num = ""
        std_inchi_key = ""
        cas_registry = ""
        
        id_elem = LightXML.find_element(compound_elem, "RegNum")
        if id_elem !== nothing
            reg_num = LightXML.content(id_elem)
            compound_data["reg_num"] = reg_num
        end
        
        # Try to find standard InChI key
        for id_elem in LightXML.get_elements_by_tagname(compound_elem, "sStandardInChI")
            std_inchi_key = LightXML.content(id_elem)
            compound_data["inchi"] = std_inchi_key
            break
        end
        
        # Try to find CAS registry number
        for id_elem in LightXML.get_elements_by_tagname(compound_elem, "CASRegistryNumber")
            cas_registry = LightXML.content(id_elem)
            compound_data["cas"] = cas_registry
            break
        end
        
        # Extract compound name
        name_elem = LightXML.find_element(compound_elem, "StandardName")
        species_name = ""
        
        if name_elem !== nothing
            species_name = LightXML.content(name_elem)
        else
            name_elem = LightXML.find_element(compound_elem, "Name")
            if name_elem !== nothing
                species_name = LightXML.content(name_elem)
            end
        end
        
        compound_data["species_name"] = species_name
        
        # Extract chemical formula
        formula_elem = LightXML.find_element(compound_elem, "ChemicalFormula")
        formula = ""
        
        if formula_elem !== nothing
            formula = LightXML.content(formula_elem)
        end
        
        compound_data["formula"] = formula
        
        # Extract molecular weight
        mw = 0.0
        mw_elem = LightXML.find_element(compound_elem, "MolecularWeight")
        
        if mw_elem !== nothing
            val_elem = LightXML.find_element(mw_elem, "nValue")
            if val_elem !== nothing
                try
                    mw = parse(Float64, LightXML.content(val_elem))
                catch e
                    @warn "Failed to parse molecular weight for $species_name: $e"
                end
            end
        end
        
        compound_data["molecular_weight"] = mw
        
        # Use compound ID or name as key
        key = isempty(reg_num) ? species_name : reg_num
        
        if !isempty(key)
            compound_data["metadata"] = Dict(
                "cas" => cas_registry,
                "inchi" => std_inchi_key,
                "source" => "ThermoML",
                "citation" => citation,
                "version" => version
            )
            
            species_data[key] = compound_data
        end
    end
    
    # Extract property data
    for prop_elem in LightXML.get_elements_by_tagname(root, "Property")
        # Get property metadata
        prop_name = ""
        prop_method = ""
        
        # Property name
        name_elem = LightXML.find_element(prop_elem, "Property-MethodID")
        if name_elem !== nothing
            prop_name_elem = LightXML.find_element(name_elem, "PropertyName")
            if prop_name_elem !== nothing
                prop_name = LightXML.content(prop_name_elem)
            end
            
            method_elem = LightXML.find_element(name_elem, "MethodName")
            if method_elem !== nothing
                prop_method = LightXML.content(method_elem)
            end
        end
        
        # Skip if not a thermodynamic property we're interested in
        if !(occursin("heat capacity", lowercase(prop_name)) || 
             occursin("enthalpy", lowercase(prop_name)) || 
             occursin("entropy", lowercase(prop_name)) || 
             occursin("gibbs", lowercase(prop_name)))
            continue
        end
        
        # Extract property data points
        for point_elem in LightXML.get_elements_by_tagname(prop_elem, "PropData")
            compound_idx = 0
            
            # Find which compound this data is for
            idx_elem = LightXML.find_element(point_elem, "nPropNumber")
            if idx_elem !== nothing
                try
                    compound_idx = parse(Int, LightXML.content(idx_elem))
                catch e
                    @warn "Failed to parse compound index: $e"
                    continue
                end
            end
            
            # Skip if no compound index found
            if compound_idx == 0
                continue
            end
            
            # Get compound ID
            compound_id = ""
            reg_num_elem = LightXML.find_element(point_elem, "RegNum")
            if reg_num_elem !== nothing
                compound_id = LightXML.content(reg_num_elem)
            end
            
            # Find compound in our data
            compound_data = nothing
            if !isempty(compound_id) && haskey(species_data, compound_id)
                compound_data = species_data[compound_id]
            else
                # Try to find by index (less reliable)
                for (_, data) in species_data
                    if get(data, "index", 0) == compound_idx
                        compound_data = data
                        break
                    end
                end
            end
            
            # Skip if compound not found
            if compound_data === nothing
                continue
            end
            
            # Initialize property data if needed
            if !haskey(compound_data, "properties")
                compound_data["properties"] = Dict{String, Any}()
            end
            
            # Initialize property type if needed
            if !haskey(compound_data["properties"], prop_name)
                compound_data["properties"][prop_name] = Dict{String, Any}(
                    "method" => prop_method,
                    "data_points" => []
                )
            end
            
            # Extract data point values
            data_point = Dict{String, Any}()
            
            # Temperature
            temp = 0.0
            temp_elem = LightXML.find_element(point_elem, "Temperature")
            if temp_elem !== nothing
                val_elem = LightXML.find_element(temp_elem, "nValue")
                if val_elem !== nothing
                    try
                        temp = parse(Float64, LightXML.content(val_elem))
                    catch e
                        @warn "Failed to parse temperature: $e"
                        continue
                    end
                end
            end
            
            data_point["temperature"] = temp
            
            # Skip if no temperature data
            if temp == 0.0
                continue
            end
            
            # Property value
            value = 0.0
            uncertainty = 0.0
            
            val_elem = LightXML.find_element(point_elem, "nValue")
            if val_elem !== nothing
                try
                    value = parse(Float64, LightXML.content(val_elem))
                catch e
                    @warn "Failed to parse property value: $e"
                    continue
                end
            end
            
            data_point["value"] = value
            
            # Uncertainty if available
            uncert_elem = LightXML.find_element(point_elem, "nUncertainty")
            if uncert_elem !== nothing
                try
                    uncertainty = parse(Float64, LightXML.content(uncert_elem))
                    data_point["uncertainty"] = uncertainty
                catch e
                    @warn "Failed to parse uncertainty: $e"
                end
            end
            
            # Units
            unit_elem = LightXML.find_element(point_elem, "Unit")
            if unit_elem !== nothing
                data_point["units"] = LightXML.content(unit_elem)
            end
            
            # Add data point to property data
            push!(compound_data["properties"][prop_name]["data_points"], data_point)
        end
    end
    
    # Process property data into tabular form
    for (key, compound_data) in species_data
        if haskey(compound_data, "properties")
            tabular_data = Dict{String, Vector{Float64}}(
                "T" => Float64[]
            )
            
            for (prop_name, prop_data) in compound_data["properties"]
                # Map ThermoML property names to our standard names
                prop_key = ""
                if occursin("heat capacity", lowercase(prop_name))
                    prop_key = "Cp"
                elseif occursin("enthalpy", lowercase(prop_name))
                    prop_key = "H"
                elseif occursin("entropy", lowercase(prop_name))
                    prop_key = "S"
                elseif occursin("gibbs", lowercase(prop_name))
                    prop_key = "G"
                else
                    continue
                end
                
                # Initialize property vectors
                tabular_data[prop_key] = Float64[]
                tabular_data["$(prop_key)_uncertainty"] = Float64[]
                
                # Add data points
                for point in prop_data["data_points"]
                    # Add temperature if not already in the list
                    temp = point["temperature"]
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
                    
                    tabular_data[prop_key][idx] = point["value"]
                    
                    # Add uncertainty if available
                    if haskey(point, "uncertainty")
                        if length(tabular_data["$(prop_key)_uncertainty"]) < idx
                            append!(tabular_data["$(prop_key)_uncertainty"], 
                                   fill(NaN, idx - length(tabular_data["$(prop_key)_uncertainty"])))
                        end
                        
                        tabular_data["$(prop_key)_uncertainty"][idx] = point["uncertainty"]
                    end
                end
            end
            
            # Get temperature range
            temps = tabular_data["T"]
            temp_min = minimum(temps)
            temp_max = maximum(temps)
            
            # Add tabular data to compound data
            compound_data["tabular_data"] = tabular_data
            compound_data["temperature_range"] = [temp_min, temp_max]
        end
    end
    
    return species_data
end

"""
    thermoml_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a ThermoML file and store the data in the database.
"""
function thermoml_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_thermoml_file(file_path)
    
    count = 0
    
    # Transaction for bulk inserts
    transaction(function()
        for (_, data) in species_data
            # Skip if no tabular data
            if !haskey(data, "tabular_data") || isempty(data["tabular_data"])
                continue
            end
            
            # Get or create species
            species_name = data["species_name"]
            formula = data["formula"]
            cas = get(get(data, "metadata", Dict()), "cas", "")
            mw = data["molecular_weight"]
            
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