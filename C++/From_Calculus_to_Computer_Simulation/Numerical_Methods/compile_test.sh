#!/bin/bash

# Compile test.cpp with all numerical methods source files
# Includes: root finding, linear systems, interpolation, ODE solvers, and optimization

# Compiler flags
CXX_FLAGS="-std=c++17 -Wall -Wextra -O2 -I."

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

# Interpolation source files
INTERPOLATION_SOURCES=(
    "src/interpolation/linear_interpolation.cpp"
    "src/interpolation/lagrange_polynomial.cpp"
    "src/interpolation/newton_divided_difference.cpp"
    "src/interpolation/cubic_spline.cpp"
)

# ODE solver source files
ODE_SOURCES=(
    "src/ode_solver/euler.cpp"
    "src/ode_solver/rk4.cpp"
    "src/ode_solver/rk45.cpp"
    "src/ode_solver/rk78.cpp"
    "src/ode_solver/bdf.cpp"
    "src/ode_solver/radau.cpp"
)

# Optimization source files
OPTIMIZATION_SOURCES=(
    "src/optimization/steepest_descent.cpp"
    "src/optimization/bfgs.cpp"
    "src/optimization/simplex.cpp"
    "src/optimization/sqp.cpp"
    "src/optimization/lagrange_multipliers.cpp"
    "src/optimization/genetic_algorithm.cpp"
    "src/optimization/pso.cpp"
)

# Combine all source files into an array
ALL_SOURCES=("test.cpp")
for file in "${ROOT_SOURCES[@]}" "${LINEAR_SOURCES[@]}" "${INTERPOLATION_SOURCES[@]}" "${ODE_SOURCES[@]}" "${OPTIMIZATION_SOURCES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Warning: Source file not found: $file"
    else
        ALL_SOURCES+=("$file")
    fi
done

# Compile
echo "Compiling test.cpp with all numerical methods..."
echo "Using compiler flags: $CXX_FLAGS"
echo ""
echo "Source files breakdown:"
echo "  - Root finding: ${#ROOT_SOURCES[@]} files"
echo "  - Linear systems: ${#LINEAR_SOURCES[@]} files"
echo "  - Interpolation: ${#INTERPOLATION_SOURCES[@]} files"
echo "  - ODE solvers: ${#ODE_SOURCES[@]} files"
echo "  - Optimization: ${#OPTIMIZATION_SOURCES[@]} files"
echo "  - Total: ${#ALL_SOURCES[@]} files (including test.cpp)"
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
