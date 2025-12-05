import numpy as np

# Perform Mann-Whitney U tests and collect p-values
p_values = [0.005497656356,
0.3030769369,
0.8342216787,
0.5939160969,
0.3492694525,
0.04392542782,
0.0004261974179,
0.00005000539841,
0.01825779853,
0.1100227956,
0.03242960796,
0.7834033132,
0.9263802457,
0.0754501538,
0.015531649]

# Bonferroni correction
alpha = 0.05
corrected_alpha = alpha / 15

# Apply the correction and print results
corrected_p_values = np.array(p_values) * 15  # Bonferroni adjustment
corrected_p_values = np.minimum(corrected_p_values, 1)  # Cap at 1.0

print("Original p-values:", p_values)
print("Corrected p-values:", corrected_p_values)
print("Significance threshold after Bonferroni correction:", corrected_alpha)

# Determine which tests are significant after correction
significant_tests = corrected_p_values < alpha
print("Significant tests after Bonferroni correction:", significant_tests)