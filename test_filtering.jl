#!/usr/bin/env julia

# Test script to verify source filtering

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData

println("Testing source filtering for NO+")
println("================================")

# Get all sources (including all theoretical)
println("\nAll sources (include_all_theoretical=true):")
all_sources = JThermodynamicsData.list_available_sources("NO+", true)
for (src_name, priority, reliability) in all_sources
    println("  - $(src_name) (priority: $(priority), reliability: $(reliability))")
end

# Get filtered sources (only highest priority theoretical)
println("\nFiltered sources (include_all_theoretical=false):")
filtered_sources = JThermodynamicsData.list_available_sources("NO+", false)
for (src_name, priority, reliability) in filtered_sources
    println("  - $(src_name) (priority: $(priority), reliability: $(reliability))")
end

println("\nTesting function implementation in json_storage.jl")
println("=================================================")

# Use the same improved filtering logic as in json_storage.jl
theoretical_sources = filter(s -> s[2] <= 4 || 
                               startswith(lowercase(s[1]), "theoretical") || 
                               lowercase(s[1]) == "theoretical" ||
                               lowercase(s[1]) == "quantum" ||
                               lowercase(s[1]) == "statistical" ||
                               lowercase(s[1]) == "quantum-statistical" ||
                               lowercase(s[1]) == "stat-thermo" ||
                               lowercase(s[1]) == "benson-group" ||
                               lowercase(s[1]) == "group-contribution", 
                        all_sources)
println("\nTheoretical sources:")
for src in theoretical_sources
    println("  - $(src[1]) (priority: $(src[2]))")
end

experimental_sources = filter(s -> s[2] > 4 && 
                               !startswith(lowercase(s[1]), "theoretical") && 
                               lowercase(s[1]) != "theoretical" &&
                               lowercase(s[1]) != "quantum" &&
                               lowercase(s[1]) != "statistical" &&
                               lowercase(s[1]) != "quantum-statistical" &&
                               lowercase(s[1]) != "stat-thermo" &&
                               lowercase(s[1]) != "benson-group" &&
                               lowercase(s[1]) != "group-contribution",
                        all_sources)
println("\nExperimental sources:")
for src in experimental_sources
    println("  - $(src[1]) (priority: $(src[2]))")
end

# Get only the highest priority theoretical source
best_theoretical = isempty(theoretical_sources) ? [] : [theoretical_sources[1]]
println("\nBest theoretical source:")
for src in best_theoretical
    println("  - $(src[1]) (priority: $(src[2]))")
end

# Combine highest priority theoretical with all experimental sources
manual_filtered_sources = vcat(best_theoretical, experimental_sources)
# Sort again to ensure correct order
sort!(manual_filtered_sources, by=s->s[2], rev=true)

println("\nManually filtered sources:")
for (src_name, priority, reliability) in manual_filtered_sources
    println("  - $(src_name) (priority: $(priority), reliability: $(reliability))")
end