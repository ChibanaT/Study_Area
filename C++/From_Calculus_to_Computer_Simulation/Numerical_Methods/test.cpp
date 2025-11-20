#include "include/root.hpp"
#include "include/linear_system.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

using namespace numerical;

// Test functions for root finding
double f1(double x) { return x * x - 4.0; }  // Root at x = 2.0
double df1(double x) { return 2.0 * x; }
double d2f1(double x) { return 2.0; }

double f2(double x) { return x * x * x - 2.0 * x - 5.0; }  // Root near x ≈ 2.0946
double df2(double x) { return 3.0 * x * x - 2.0; }
double d2f2(double x) { return 6.0 * x; }

double f3(double x) { return std::exp(x) - 3.0 * x; }  // Root near x ≈ 0.619
double df3(double x) { return std::exp(x) - 3.0; }
double d2f3(double x) { return std::exp(x); }

// Helper function to print root result
void print_root_result(const std::string& method_name, const RootResult& result, double expected_root, double tolerance = 1e-6) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Root: " << result.root << " (expected: " << expected_root << ")\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    std::cout << "  Error: " << std::fabs(result.root - expected_root) << "\n";
    if (std::fabs(result.root - expected_root) < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

// Helper function to print linear system result
void print_linear_result(const std::string& method_name, const LinearResult& result, 
                         const std::vector<double>& expected, double tolerance = 1e-6) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Solution: [";
    for (size_t i = 0; i < result.solution.size(); ++i) {
        std::cout << result.solution[i];
        if (i < result.solution.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
    std::cout << "  Expected: [";
    for (size_t i = 0; i < expected.size(); ++i) {
        std::cout << expected[i];
        if (i < expected.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    
    double max_error = 0.0;
    for (size_t i = 0; i < result.solution.size() && i < expected.size(); ++i) {
        double error = std::fabs(result.solution[i] - expected[i]);
        if (error > max_error) max_error = error;
    }
    std::cout << "  Max Error: " << max_error << "\n";
    if (max_error < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

// Test all root finding methods
void test_all_root_methods() {
    std::cout << "\n========================================\n";
    std::cout << "TESTING ALL ROOT FINDING METHODS\n";
    std::cout << "========================================\n";

    double tol = 1e-10;
    int max_iter = 1000;

    // Test 1: f(x) = x^2 - 4, root at x = 2.0
    std::cout << "\n--- Test 1: f(x) = x^2 - 4 (root at x = 2.0) ---\n";
    
    RootResult result;
    
    try {
        result = bisect(f1, 0.0, 5.0, tol, max_iter);
        print_root_result("Bisection", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Bisection Error: " << e.what() << "\n";
    }
    
    try {
        result = false_position(f1, 0.0, 5.0, tol, max_iter);
        print_root_result("False Position", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_false_position(f1, 0.0, 5.0, tol, max_iter);
        print_root_result("Modified False Position", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Modified False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = secant(f1, 0.0, 5.0, tol, max_iter);
        print_root_result("Secant", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Secant Error: " << e.what() << "\n";
    }
    
    try {
        result = newton_raphson(f1, df1, 3.0, tol, max_iter);
        print_root_result("Newton-Raphson", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Newton-Raphson Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_newton(f1, df1, 3.0, tol, max_iter);
        print_root_result("Modified Newton", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Modified Newton Error: " << e.what() << "\n";
    }
    
    try {
        result = halley(f1, df1, d2f1, 3.0, tol, max_iter);
        print_root_result("Halley", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Halley Error: " << e.what() << "\n";
    }
    
    try {
        result = brent(f1, 0.0, 5.0, tol, max_iter);
        print_root_result("Brent", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Brent Error: " << e.what() << "\n";
    }
    
    try {
        result = muller(f1, 0.0, 2.5, 5.0, tol, max_iter);
        print_root_result("Muller", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Muller Error: " << e.what() << "\n";
    }
    
    // Bairstow test (polynomial: x^2 - 4 = 0, coefficients: [-4, 0, 1])
    try {
        std::vector<double> coeffs = {-4.0, 0.0, 1.0};
        result = bairstow(coeffs, 1.0, 1.0, tol, max_iter);
        print_root_result("Bairstow", result, 2.0);
    } catch (const std::exception& e) {
        std::cout << "  Bairstow Error: " << e.what() << "\n";
    }

    // Test 2: f(x) = x^3 - 2x - 5, root near x ≈ 2.0946
    std::cout << "\n--- Test 2: f(x) = x^3 - 2x - 5 (root near x ≈ 2.0946) ---\n";
    
    RootResult result;
    
    try {
        result = bisect(f2, 0.0, 5.0, tol, max_iter);
        print_root_result("Bisection", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Bisection Error: " << e.what() << "\n";
    }
    
    try {
        result = false_position(f2, 0.0, 5.0, tol, max_iter);
        print_root_result("False Position", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_false_position(f2, 0.0, 5.0, tol, max_iter);
        print_root_result("Modified False Position", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Modified False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = secant(f2, 0.0, 5.0, tol, max_iter);
        print_root_result("Secant", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Secant Error: " << e.what() << "\n";
    }
    
    try {
        result = newton_raphson(f2, df2, 3.0, tol, max_iter);
        print_root_result("Newton-Raphson", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Newton-Raphson Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_newton(f2, df2, 3.0, tol, max_iter);
        print_root_result("Modified Newton", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Modified Newton Error: " << e.what() << "\n";
    }
    
    try {
        result = halley(f2, df2, d2f2, 3.0, tol, max_iter);
        print_root_result("Halley", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Halley Error: " << e.what() << "\n";
    }
    
    try {
        result = brent(f2, 0.0, 5.0, tol, max_iter);
        print_root_result("Brent", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Brent Error: " << e.what() << "\n";
    }
    
    try {
        result = muller(f2, 0.0, 2.5, 5.0, tol, max_iter);
        print_root_result("Muller", result, 2.0945514815);
    } catch (const std::exception& e) {
        std::cout << "  Muller Error: " << e.what() << "\n";
    }

    // Test 3: f(x) = exp(x) - 3x, root near x ≈ 0.619
    std::cout << "\n--- Test 3: f(x) = exp(x) - 3x (root near x ≈ 0.619) ---\n";
    
    RootResult result;
    
    try {
        result = bisect(f3, -1.0, 2.0, tol, max_iter);
        print_root_result("Bisection", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Bisection Error: " << e.what() << "\n";
    }
    
    try {
        result = false_position(f3, -1.0, 2.0, tol, max_iter);
        print_root_result("False Position", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_false_position(f3, -1.0, 2.0, tol, max_iter);
        print_root_result("Modified False Position", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Modified False Position Error: " << e.what() << "\n";
    }
    
    try {
        result = secant(f3, -1.0, 2.0, tol, max_iter);
        print_root_result("Secant", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Secant Error: " << e.what() << "\n";
    }
    
    try {
        result = newton_raphson(f3, df3, 1.0, tol, max_iter);
        print_root_result("Newton-Raphson", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Newton-Raphson Error: " << e.what() << "\n";
    }
    
    try {
        result = modified_newton(f3, df3, 1.0, tol, max_iter);
        print_root_result("Modified Newton", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Modified Newton Error: " << e.what() << "\n";
    }
    
    try {
        result = halley(f3, df3, d2f3, 1.0, tol, max_iter);
        print_root_result("Halley", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Halley Error: " << e.what() << "\n";
    }
    
    try {
        result = brent(f3, -1.0, 2.0, tol, max_iter);
        print_root_result("Brent", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Brent Error: " << e.what() << "\n";
    }
    
    try {
        result = muller(f3, -1.0, 0.5, 2.0, tol, max_iter);
        print_root_result("Muller", result, 0.6190612867);
    } catch (const std::exception& e) {
        std::cout << "  Muller Error: " << e.what() << "\n";
    }
}

// Test all linear system methods
void test_all_linear_system_methods() {
    std::cout << "\n\n========================================\n";
    std::cout << "TESTING ALL LINEAR SYSTEM METHODS\n";
    std::cout << "========================================\n";

    double tol = 1e-10;
    int max_iter = 1000;

    // Test 1: Simple 2x2 system
    // 2x + y = 5
    // x + 3y = 10
    // Solution: x = 1, y = 3
    std::cout << "\n--- Test 1: 2x2 System (x=1, y=3) ---\n";
    
    std::vector<std::vector<double>> A1 = {
        {2.0, 1.0},
        {1.0, 3.0}
    };
    std::vector<double> b1 = {5.0, 10.0};
    std::vector<double> expected1 = {1.0, 3.0};

    LinearResult result;
    
    try {
        result = cramer(A1, b1);
        print_linear_result("Cramer", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Cramer Error: " << e.what() << "\n";
    }
    
    try {
        result = gaussian_elimination(A1, b1);
        print_linear_result("Gaussian Elimination", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Gaussian Elimination Error: " << e.what() << "\n";
    }
    
    try {
        result = pivot_gauss(A1, b1);
        print_linear_result("Pivoted Gaussian", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Pivoted Gaussian Error: " << e.what() << "\n";
    }
    
    try {
        result = gauss_jordan(A1, b1);
        print_linear_result("Gauss-Jordan", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Jordan Error: " << e.what() << "\n";
    }
    
    try {
        result = lu_decomposition(A1, b1);
        print_linear_result("LU Decomposition", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  LU Decomposition Error: " << e.what() << "\n";
    }
    
    // Cholesky requires symmetric positive definite matrix - skip for A1
    std::cout << "\n  Cholesky Decomposition: Skipped (matrix not symmetric)\n";
    
    try {
        result = jacobi(A1, b1, max_iter, tol);
        print_linear_result("Jacobi", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Jacobi Error: " << e.what() << "\n";
    }
    
    try {
        result = gauss_seidel(A1, b1, max_iter, tol);
        print_linear_result("Gauss-Seidel", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Seidel Error: " << e.what() << "\n";
    }
    
    try {
        result = sor(A1, b1, max_iter, tol, 1.2);
        print_linear_result("SOR (omega=1.2)", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  SOR Error: " << e.what() << "\n";
    }
    
    // Conjugate Gradient requires symmetric positive definite matrix - skip for A1
    std::cout << "\n  Conjugate Gradient: Skipped (matrix not symmetric)\n";

    // Test 2: 3x3 system
    // x + 2y + 3z = 14
    // 2x + 5y + 2z = 18
    // 3x + y + 5z = 20
    // Solution: x = 1, y = 2, z = 3
    std::cout << "\n--- Test 2: 3x3 System (x=1, y=2, z=3) ---\n";
    
    std::vector<std::vector<double>> A2 = {
        {1.0, 2.0, 3.0},
        {2.0, 5.0, 2.0},
        {3.0, 1.0, 5.0}
    };
    std::vector<double> b2 = {14.0, 18.0, 20.0};
    std::vector<double> expected2 = {1.0, 2.0, 3.0};

    LinearResult result;
    
    try {
        result = gaussian_elimination(A2, b2);
        print_linear_result("Gaussian Elimination", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Gaussian Elimination Error: " << e.what() << "\n";
    }
    
    try {
        result = pivot_gauss(A2, b2);
        print_linear_result("Pivoted Gaussian", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Pivoted Gaussian Error: " << e.what() << "\n";
    }
    
    try {
        result = gauss_jordan(A2, b2);
        print_linear_result("Gauss-Jordan", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Jordan Error: " << e.what() << "\n";
    }
    
    try {
        result = lu_decomposition(A2, b2);
        print_linear_result("LU Decomposition", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  LU Decomposition Error: " << e.what() << "\n";
    }
    
    try {
        result = jacobi(A2, b2, max_iter, tol);
        print_linear_result("Jacobi", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Jacobi Error: " << e.what() << "\n";
    }
    
    try {
        result = gauss_seidel(A2, b2, max_iter, tol);
        print_linear_result("Gauss-Seidel", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Seidel Error: " << e.what() << "\n";
    }
    
    try {
        result = sor(A2, b2, max_iter, tol, 1.1);
        print_linear_result("SOR (omega=1.1)", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  SOR Error: " << e.what() << "\n";
    }
    
    // Conjugate Gradient requires symmetric positive definite matrix - skip for A2
    std::cout << "\n  Conjugate Gradient: Skipped (matrix not symmetric)\n";

    // Test 3: Tridiagonal system for Thomas algorithm
    // 2x + y = 3
    // x + 2y + z = 6
    // y + 2z = 7
    // Solution: x = 1, y = 1, z = 3
    std::cout << "\n--- Test 3: Tridiagonal System (x=1, y=1, z=3) ---\n";
    
    std::vector<std::vector<double>> A3 = {
        {2.0, 1.0, 0.0},
        {1.0, 2.0, 1.0},
        {0.0, 1.0, 2.0}
    };
    std::vector<double> b3 = {3.0, 6.0, 7.0};
    std::vector<double> expected3 = {1.0, 1.0, 3.0};

    LinearResult result;
    
    try {
        result = thomas_algorithm(A3, b3);
        print_linear_result("Thomas Algorithm", result, expected3);
    } catch (const std::exception& e) {
        std::cout << "  Thomas Algorithm Error: " << e.what() << "\n";
    }
    
    try {
        result = gaussian_elimination(A3, b3);
        print_linear_result("Gaussian Elimination", result, expected3);
    } catch (const std::exception& e) {
        std::cout << "  Gaussian Elimination Error: " << e.what() << "\n";
    }
    
    try {
        result = lu_decomposition(A3, b3);
        print_linear_result("LU Decomposition", result, expected3);
    } catch (const std::exception& e) {
        std::cout << "  LU Decomposition Error: " << e.what() << "\n";
    }
    
    // A3 is symmetric, test Cholesky and Conjugate Gradient
    try {
        result = cholesky_decomposition(A3, b3);
        print_linear_result("Cholesky Decomposition", result, expected3);
    } catch (const std::exception& e) {
        std::cout << "  Cholesky Decomposition Error: " << e.what() << "\n";
    }
    
    try {
        result = conjugate_gradient(A3, b3, max_iter, tol);
        print_linear_result("Conjugate Gradient", result, expected3);
    } catch (const std::exception& e) {
        std::cout << "  Conjugate Gradient Error: " << e.what() << "\n";
    }

    // Test 4: Symmetric positive definite matrix for Cholesky
    // 4x + 2y = 10
    // 2x + 5y = 16
    // Solution: x = 1.125, y = 2.75
    std::cout << "\n--- Test 4: SPD System for Cholesky (x=1.125, y=2.75) ---\n";
    
    std::vector<std::vector<double>> A4 = {
        {4.0, 2.0},
        {2.0, 5.0}
    };
    std::vector<double> b4 = {10.0, 16.0};
    std::vector<double> expected4 = {1.125, 2.75};

    LinearResult result;
    
    try {
        result = cholesky_decomposition(A4, b4);
        print_linear_result("Cholesky Decomposition", result, expected4);
    } catch (const std::exception& e) {
        std::cout << "  Cholesky Decomposition Error: " << e.what() << "\n";
    }
    
    try {
        result = conjugate_gradient(A4, b4, max_iter, tol);
        print_linear_result("Conjugate Gradient", result, expected4);
    } catch (const std::exception& e) {
        std::cout << "  Conjugate Gradient Error: " << e.what() << "\n";
    }
    
    // Also test other methods on this SPD matrix
    try {
        result = gaussian_elimination(A4, b4);
        print_linear_result("Gaussian Elimination", result, expected4);
    } catch (const std::exception& e) {
        std::cout << "  Gaussian Elimination Error: " << e.what() << "\n";
    }
    
    try {
        result = lu_decomposition(A4, b4);
        print_linear_result("LU Decomposition", result, expected4);
    } catch (const std::exception& e) {
        std::cout << "  LU Decomposition Error: " << e.what() << "\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "NUMERICAL METHODS TEST SUITE\n";
    std::cout << "========================================\n";

    // Test all root finding methods
    test_all_root_methods();

    // Test all linear system methods
    test_all_linear_system_methods();

    std::cout << "\n\n========================================\n";
    std::cout << "TESTING COMPLETE\n";
    std::cout << "========================================\n";

    return 0;
}

