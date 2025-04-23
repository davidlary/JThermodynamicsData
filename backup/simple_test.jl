# Simple test script for specific functionality

function is_ionic_species_check(species_name, formula)
    return endswith(species_name, "+") || endswith(species_name, "-") ||
           endswith(formula, "+") || endswith(formula, "-") ||
           species_name == "e-" || formula == "e-"
end

# Test our changes on example species
species_tests = [
    ("He", "He", false),
    ("He+", "He+", true),
    ("Xe", "Xe", false),
    ("Xe+", "Xe+", true),
    ("N2", "N2", false),
    ("NO+", "NO+", true),
    ("e-", "e-", true)
]

println("Testing ionic species detection:")
for (species, formula, expected) in species_tests
    result = is_ionic_species_check(species, formula)
    status = result == expected ? "✓" : "✗"
    println("  $status $species, $formula => is_ionic: $result (expected: $expected)")
end

# Check the plots directory to see what ion plots are available
println("\nFound ionic species plots:")
ion_plots = filter(f -> (endswith(f, ".png") && (contains(f, "+") || contains(f, "-"))), 
                  readdir(joinpath(@__DIR__, "plots")))

for plot in ion_plots
    println("  $plot")
end