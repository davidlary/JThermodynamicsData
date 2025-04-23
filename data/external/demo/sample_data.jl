#!/usr/bin/env julia

# Sample data for hierarchical thermodynamic calculator demonstration

"""
This file contains sample data for the hierarchical thermodynamic calculator demonstration.
It includes NASA-7 polynomial coefficients for common species from multiple sources
for testing the hierarchical loading and prioritization of data.
"""

# Dictionary structure: source -> species -> properties
SAMPLE_DATA = Dict(
    "theoretical" => Dict(
        "N2" => Dict(
            "formula" => "N2",
            "name" => "Nitrogen",
            "cas" => "7727-37-9",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.53100, -0.000123, -5.02999e-7, 2.43530e-9, -1.40881e-12, -1046.97, 2.96747],  # 200-1000K
                [2.95257, 0.0013969, -4.92631e-7, 7.86010e-11, -4.60755e-15, -923.93, 5.87189]   # 1000-6000K
            ],
            "reliability" => 2.5,
            "priority" => 1
        ),
        "O2" => Dict(
            "formula" => "O2",
            "name" => "Oxygen",
            "cas" => "7782-44-7",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.78246, -0.00299673, 9.84730e-6, -9.68129e-9, 3.24372e-12, -1063.94, 3.65768],  # 200-1000K
                [3.66096, 0.000656366, -1.41149e-7, 2.05798e-11, -1.29913e-15, -1215.98, 3.41536] # 1000-6000K
            ],
            "reliability" => 2.5,
            "priority" => 1
        ),
        "H2O" => Dict(
            "formula" => "H2O",
            "name" => "Water",
            "cas" => "7732-18-5",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [4.19864, -0.00203643, 6.52040e-6, -5.48797e-9, 1.77197e-12, -30293.7, -0.849032], # 200-1000K
                [2.67704, 0.00297318, -7.73769e-7, 9.44335e-11, -4.26900e-15, -29885.9, 6.88255]  # 1000-6000K
            ],
            "reliability" => 2.5,
            "priority" => 1
        )
    ),
    "gri-mech" => Dict(
        "N2" => Dict(
            "formula" => "N2",
            "name" => "Nitrogen",
            "cas" => "7727-37-9",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.29867, 0.00140824, -3.96322e-6, 5.64152e-9, -2.44485e-12, -1020.9, 3.95037],  # 200-1000K
                [2.92664, 0.00148798, -5.68476e-7, 1.00970e-10, -6.75335e-15, -922.798, 5.98053] # 1000-6000K
            ],
            "reliability" => 3.5,
            "priority" => 3
        ),
        "CH4" => Dict(
            "formula" => "CH4",
            "name" => "Methane",
            "cas" => "74-82-8",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [5.14987, -0.0136709, 4.91800e-5, -4.84743e-8, 1.66693e-11, -10246.6, -4.64130],  # 200-1000K
                [1.65326, 0.0100263, -3.31661e-6, 5.36483e-10, -3.14696e-14, -10009.6, 9.99486]   # 1000-6000K
            ],
            "reliability" => 3.5,
            "priority" => 3
        )
    ),
    "burcat" => Dict(
        "N2" => Dict(
            "formula" => "N2",
            "name" => "Nitrogen",
            "cas" => "7727-37-9",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.53100, -0.000123, -5.02999e-7, 2.43530e-9, -1.40881e-12, -1046.97, 2.96747],  # 200-1000K
                [2.95257, 0.0013969, -4.92631e-7, 7.86010e-11, -4.60755e-15, -923.93, 5.87189]   # 1000-6000K
            ],
            "reliability" => 4.2,
            "priority" => 7
        ),
        "O2" => Dict(
            "formula" => "O2",
            "name" => "Oxygen",
            "cas" => "7782-44-7",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.78246, -0.00299673, 9.84730e-6, -9.68129e-9, 3.24372e-12, -1063.94, 3.65768],  # 200-1000K
                [3.66096, 0.000656366, -1.41149e-7, 2.05798e-11, -1.29913e-15, -1215.98, 3.41536] # 1000-6000K
            ],
            "reliability" => 4.2,
            "priority" => 7
        ),
        "H2O" => Dict(
            "formula" => "H2O",
            "name" => "Water",
            "cas" => "7732-18-5",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [4.19864, -0.00203643, 6.52040e-6, -5.48797e-9, 1.77197e-12, -30293.7, -0.849032], # 200-1000K
                [2.67704, 0.00297318, -7.73769e-7, 9.44335e-11, -4.26900e-15, -29885.9, 6.88255]  # 1000-6000K
            ],
            "reliability" => 4.2,
            "priority" => 7
        )
    ),
    "janaf" => Dict(
        "N2" => Dict(
            "formula" => "N2",
            "name" => "Nitrogen",
            "cas" => "7727-37-9",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.50342, 0.000499197, -2.30471e-7, 1.04436e-9, -5.26607e-13, -1043.76, 2.96747], # 200-1000K
                [2.95257, 0.0013969, -4.92631e-7, 7.86010e-11, -4.60755e-15, -923.93, 5.87189]   # 1000-6000K
            ],
            "reliability" => 4.5,
            "priority" => 8
        ),
        "O2" => Dict(
            "formula" => "O2",
            "name" => "Oxygen",
            "cas" => "7782-44-7",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.78246, -0.00299673, 9.84730e-6, -9.68129e-9, 3.24372e-12, -1063.94, 3.65768],  # 200-1000K
                [3.66096, 0.000656366, -1.41149e-7, 2.05798e-11, -1.29913e-15, -1215.98, 3.41536] # 1000-6000K
            ],
            "reliability" => 4.5,
            "priority" => 8
        ),
        "H2O" => Dict(
            "formula" => "H2O",
            "name" => "Water",
            "cas" => "7732-18-5",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [4.19864, -0.00203643, 6.52040e-6, -5.48797e-9, 1.77197e-12, -30293.7, -0.849032], # 200-1000K
                [2.67704, 0.00297318, -7.73769e-7, 9.44335e-11, -4.26900e-15, -29885.9, 6.88255]  # 1000-6000K
            ],
            "reliability" => 4.5,
            "priority" => 8
        )
    ),
    "atct" => Dict(
        "N2" => Dict(
            "formula" => "N2",
            "name" => "Nitrogen",
            "cas" => "7727-37-9", 
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.52525, 0.000287533, -7.69284e-7, 1.07951e-9, -5.08524e-13, -1043.65, 2.96747], # 200-1000K
                [2.95257, 0.0013969, -4.92631e-7, 7.86010e-11, -4.60755e-15, -923.93, 5.87189]   # 1000-6000K
            ],
            "reliability" => 4.9,
            "priority" => 9
        ),
        "O2" => Dict(
            "formula" => "O2",
            "name" => "Oxygen",
            "cas" => "7782-44-7",
            "polynomial_type" => "nasa7",
            "temperature_ranges" => [[200.0, 1000.0], [1000.0, 6000.0]],
            "coefficients" => [
                [3.78246, -0.00299673, 9.84730e-6, -9.68129e-9, 3.24372e-12, -1063.94, 3.65768],  # 200-1000K
                [3.66096, 0.000656366, -1.41149e-7, 2.05798e-11, -1.29913e-15, -1215.98, 3.41536] # 1000-6000K
            ],
            "reliability" => 4.9,
            "priority" => 9
        )
    )
)