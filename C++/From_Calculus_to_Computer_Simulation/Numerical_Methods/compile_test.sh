#!/bin/bash

# Compile test.cpp with all root finding and linear system source files

# Compiler flags
CXX_FLAGS="-std=c++17 -Wall -Wextra -O2"

# Root finding source files
ROOT_SOURCES=(
    "src/root/bisection.cpp"
    "src/root/false_position.cpp"
    "src/root/modified_false_position.cpp"
    "src/root/secant.cpp"
    "src/root/newton_raphson.cpp"
    "src/root/modified_newton.cpp"
    "src/root/halley.cpp"
    "src/root/brent.cpp"
    "src/root/muller.cpp"
    "src/root/bairstow.cpp"
)

# Linear system source files
LINEAR_SOURCES=(
    "src/linear_system/cramer.cpp"
    "src/linear_system/gaussian_elimination.cpp"
    "src/linear_system/pivot_gauss.cpp"
    "src/linear_system/gauss_jordan.cpp"
    "src/linear_system/lu_decomposition.cpp"
    "src/linear_system/cholesky_decomposition.cpp"
    "src/linear_system/thomas_algorithm.cpp"
    "src/linear_system/jacobi.cpp"
    "src/linear_system/gauss_seidel.cpp"
    "src/linear_system/sor.cpp"
    "src/linear_system/conjugate_gradient.cpp"
)

# Combine all source files into an array
ALL_SOURCES=("test.cpp")
for file in "${ROOT_SOURCES[@]}" "${LINEAR_SOURCES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Warning: Source file not found: $file"
    else
        ALL_SOURCES+=("$file")
    fi
done

# Compile
echo "Compiling test.cpp with all source files..."
echo "Using compiler flags: $CXX_FLAGS"
echo "Total source files: ${#ALL_SOURCES[@]}"
echo ""

g++ $CXX_FLAGS "${ALL_SOURCES[@]}" -o test.exe

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Compilation successful!"
    echo "Running test.exe..."
    echo "========================================"
    echo ""
    ./test.exe
    EXIT_CODE=$?
    echo ""
    echo "========================================"
    echo "Test execution completed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
else
    echo ""
    echo "✗ Compilation failed!"
    exit 1
fi
