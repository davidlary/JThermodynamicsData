# Script to analyze plot titles in the PNG files

using Base64

function extract_image_metadata(png_file)
    # Open the PNG file
    data = open(png_file, "r") do file
        read(file)
    end
    
    # Try to extract any text that might be in the image
    # This is a very simple approach that just looks for text chunks in the PNG
    metadata = Dict()
    metadata["filename"] = basename(png_file)
    metadata["filesize"] = filesize(png_file)
    
    # Output the first 200 bytes in hex for debugging
    metadata["header"] = bytes2hex(data[1:min(200, length(data))])
    
    return metadata
end

# Define the plot files to check
plots_to_check = [
    joinpath(@__DIR__, "plots", "He+_all_properties.png"),
    joinpath(@__DIR__, "plots", "Xe+_all_properties.png"),
    joinpath(@__DIR__, "plots", "Xe_all_properties.png")
]

# Check each plot
println("Analyzing plot files to check for title information:")
for plot_file in plots_to_check
    if isfile(plot_file)
        println("\nAnalyzing: $(basename(plot_file))")
        metadata = extract_image_metadata(plot_file)
        
        # Display the metadata
        for (key, value) in metadata
            if key == "header"
                println("  Header (hex): $(value[1:60])...")
            else
                println("  $key: $value")
            end
        end
    else
        println("File not found: $plot_file")
    end
end