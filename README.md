# Elastic Hash Table

A simplified implementation of the Elastic Hashing algorithm described by Martín Farach-Colton, Andrew Krapivin and William Kuszmaul in [Optimal Bounds for Open Addressing Without Reordering" by Martín Farach-Colton, Andrew Krapivin, and William Kuszmaul (2025)](https://arxiv.org/pdf/2501.02305)

## Features

- O(1) amortized expected probe complexity
- O(log δ⁻¹) worst-case expected probe complexity
- No reordering of elements
- Configurable load factor
- Performance statistics tracking
- Thread-safe search operations
- Uses geometrically decreasing array sizes
- Implements the φ(i,j) injection function for probe sequences
- Handles batch-based insertions as described in the paper
- Maintains proper load balancing across subarrays

