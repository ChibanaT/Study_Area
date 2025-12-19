#include "include/root.hpp"
#include "include/linear_system.hpp"
#include "include/interpolation.hpp"
#include "include/ode_solver.hpp"
#include "include/optimization.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <random>

using namespace numerical;

// ============================================
// TEST FUNCTIONS FOR ROOT FINDING
// ============================================

// 1. Bisection: Ideal for functions with clear sign change
// f(x) = x^5 - x - 1, root near x ≈ 1.1673
double f_bisection(double x) { return x*x*x*x*x - x - 1.0; }

// 2. False Position: Better convergence than bisection for some cases
// f(x) = x^3 - x - 2, root near x ≈ 1.5214
double f_false_position(double x) { return x*x*x - x - 2.0; }
double df_false_position(double x) { return 3.0*x*x - 1.0; }

// 3. Modified False Position: Handles slow convergence of false position
// f(x) = x^4 - 10x^2 + 9, root at x = 3.0
double f_modified_false_position(double x) { return x*x*x*x - 10.0*x*x + 9.0; }
double df_modified_false_position(double x) { return 4.0*x*x*x - 20.0*x; }

// 4. Secant: When derivative is not available
// f(x) = x*exp(x) - 1, root near x ≈ 0.5671
double f_secant(double x) { return x * std::exp(x) - 1.0; }

// 5. Newton-Raphson: Well-behaved function with good derivative
// f(x) = x^2 - 2, root at x = sqrt(2) ≈ 1.4142
double f_newton(double x) { return x*x - 2.0; }
double df_newton(double x) { return 2.0*x; }
double d2f_newton(double x) { return 2.0; }

// 6. Modified Newton: When Newton may have convergence issues
// f(x) = x^3 - 2x - 5, root near x ≈ 2.0946
double f_modified_newton(double x) { return x*x*x - 2.0*x - 5.0; }
double df_modified_newton(double x) { return 3.0*x*x - 2.0; }
double d2f_modified_newton(double x) { return 6.0*x; }

// 7. Halley: When second derivative is available
// f(x) = exp(x) - 2, root at x = ln(2) ≈ 0.6931
double f_halley(double x) { return std::exp(x) - 2.0; }
double df_halley(double x) { return std::exp(x); }
double d2f_halley(double x) { return std::exp(x); }

// 8. Brent: Robust method for difficult problems
// f(x) = x*sin(x) - 1, root near x ≈ 1.1142
double f_brent(double x) { return x * std::sin(x) - 1.0; }

// 9. Muller: For finding multiple or complex roots
// f(x) = x^4 - 5x^2 + 4, roots at x = ±2, ±1
double f_muller(double x) { return x*x*x*x - 5.0*x*x + 4.0; }

// ============================================
// HELPER FUNCTIONS
// ============================================

void print_root_result(const std::string& method_name, const RootResult& result, double expected_root, double tolerance = 1e-6) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Root: " << result.root << " (expected: " << expected_root << ")\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    std::cout << "  Error: " << std::fabs(result.root - expected_root) << "\n";
    if (result.converged && std::fabs(result.root - expected_root) < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

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
    if (result.converged && max_error < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

// ============================================
// ROOT FINDING TESTS
// ============================================

void test_all_root_methods() {
    std::cout << "\n========================================\n";
    std::cout << "TESTING ALL ROOT FINDING METHODS\n";
    std::cout << "========================================\n";

    double tol = 1e-10;
    int max_iter = 1000;
    RootResult result;

    // Test 1: Bisection - Ideal for functions with clear sign change
    std::cout << "\n--- Test 1: Bisection - f(x) = x^5 - x - 1 (root ≈ 1.1673) ---\n";
    try {
        result = bisect(f_bisection, 1.0, 2.0, tol, max_iter);
        print_root_result("Bisection", result, 1.1673039783, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Bisection Error: " << e.what() << "\n";
    }

    // Test 2: False Position - Better convergence for some cases
    std::cout << "\n--- Test 2: False Position - f(x) = x^3 - x - 2 (root ≈ 1.5214) ---\n";
    try {
        result = false_position(f_false_position, 1.0, 2.0, tol, max_iter);
        print_root_result("False Position", result, 1.5213797068, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  False Position Error: " << e.what() << "\n";
    }

    // Test 3: Modified False Position - Handles slow convergence
    std::cout << "\n--- Test 3: Modified False Position - f(x) = x^4 - 10x^2 + 9 (root = 3.0) ---\n";
    try {
        result = modified_false_position(f_modified_false_position, 2.0, 4.0, tol, max_iter);
        print_root_result("Modified False Position", result, 3.0, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Modified False Position Error: " << e.what() << "\n";
    }

    // Test 4: Secant - When derivative is not available
    std::cout << "\n--- Test 4: Secant - f(x) = x*exp(x) - 1 (root ≈ 0.5671) ---\n";
    try {
        result = secant(f_secant, 0.0, 1.0, tol, max_iter);
        print_root_result("Secant", result, 0.5671432904, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Secant Error: " << e.what() << "\n";
    }

    // Test 5: Newton-Raphson - Well-behaved function
    std::cout << "\n--- Test 5: Newton-Raphson - f(x) = x^2 - 2 (root = sqrt(2) ≈ 1.4142) ---\n";
    try {
        result = newton_raphson(f_newton, df_newton, 1.5, tol, max_iter);
        print_root_result("Newton-Raphson", result, 1.4142135624, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Newton-Raphson Error: " << e.what() << "\n";
    }

    // Test 6: Modified Newton - When Newton may have issues
    std::cout << "\n--- Test 6: Modified Newton - f(x) = x^3 - 2x - 5 (root ≈ 2.0946) ---\n";
    try {
        result = modified_newton(f_modified_newton, df_modified_newton, 2.5, tol, max_iter);
        print_root_result("Modified Newton", result, 2.0945514815, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Modified Newton Error: " << e.what() << "\n";
    }

    // Test 7: Halley - When second derivative is available
    std::cout << "\n--- Test 7: Halley - f(x) = exp(x) - 2 (root = ln(2) ≈ 0.6931) ---\n";
    try {
        result = halley(f_halley, df_halley, d2f_halley, 0.5, tol, max_iter);
        print_root_result("Halley", result, 0.6931471806, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Halley Error: " << e.what() << "\n";
    }

    // Test 8: Brent - Robust method for difficult problems
    std::cout << "\n--- Test 8: Brent - f(x) = x*sin(x) - 1 (root ≈ 1.1142) ---\n";
    try {
        result = brent(f_brent, 0.0, 2.0, tol, max_iter);
        print_root_result("Brent", result, 1.1141571409, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Brent Error: " << e.what() << "\n";
    }

    // Test 9: Muller - For finding multiple roots
    std::cout << "\n--- Test 9: Muller - f(x) = x^4 - 5x^2 + 4 (root = 2.0) ---\n";
    try {
        result = muller(f_muller, 1.5, 1.8, 2.2, tol, max_iter);
        print_root_result("Muller", result, 2.0, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Muller Error: " << e.what() << "\n";
    }

    // Test 10: Bairstow - For polynomials (x^2 - 4 = 0, root = 2.0)
    std::cout << "\n--- Test 10: Bairstow - x^2 - 4 = 0 (root = 2.0) ---\n";
    try {
        // Coefficients in order [a0, a1, a2] for a0 + a1*x + a2*x^2
        // For x^2 - 4 = 0: a0 = -4, a1 = 0, a2 = 1
        std::vector<double> coeffs = {-4.0, 0.0, 1.0};
        result = bairstow(coeffs, 0.1, -3.9, tol, max_iter);
        // Bairstow finds quadratic factor x^2 - r*x - s = 0
        // Extract the root closest to 2.0
        if (result.converged && !result.history.empty()) {
            // result.root is one root, history[0] might be the other
            double root1 = result.root;
            double root2 = (result.history.size() > 0) ? result.history[0] : root1;
            result.root = (std::fabs(root1 - 2.0) < std::fabs(root2 - 2.0)) ? root1 : root2;
        }
        print_root_result("Bairstow", result, 2.0, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  Bairstow Error: " << e.what() << "\n";
    }
}

// ============================================
// LINEAR SYSTEM TESTS
// ============================================

void test_all_linear_system_methods() {
    std::cout << "\n\n========================================\n";
    std::cout << "TESTING ALL LINEAR SYSTEM METHODS\n";
    std::cout << "========================================\n";

    double tol = 1e-10;
    int max_iter = 1000;
    LinearResult result;

    // Test 1: Cramer - Ideal for small systems (2x2)
    std::cout << "\n--- Test 1: Cramer - Small 2x2 System (x=1, y=3) ---\n";
    std::vector<std::vector<double>> A1 = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<double> b1 = {5.0, 10.0};
    std::vector<double> expected1 = {1.0, 3.0};
    try {
        result = cramer(A1, b1);
        print_linear_result("Cramer", result, expected1);
    } catch (const std::exception& e) {
        std::cout << "  Cramer Error: " << e.what() << "\n";
    }

    // Test 2: Gaussian Elimination - Standard method
    std::cout << "\n--- Test 2: Gaussian Elimination - Standard 3x3 System (x=1, y=2, z=3) ---\n";
    std::vector<std::vector<double>> A2 = {
        {1.0, 2.0, 3.0},
        {2.0, 5.0, 2.0},
        {3.0, 1.0, 5.0}
    };
    std::vector<double> b2 = {14.0, 18.0, 20.0};
    std::vector<double> expected2 = {1.0, 2.0, 3.0};
    try {
        result = gaussian_elimination(A2, b2);
        print_linear_result("Gaussian Elimination", result, expected2);
    } catch (const std::exception& e) {
        std::cout << "  Gaussian Elimination Error: " << e.what() << "\n";
    }

    // Test 3: Pivoted Gaussian - System with small pivots (near-singular)
    std::cout << "\n--- Test 3: Pivoted Gaussian - Near-Singular System ---\n";
    // System: 0.001x + y + z = 2.001
    //         x + 2y + z = 4
    //         2x + y + 3z = 6
    // Solution: x = 1, y = 1, z = 1
    std::vector<std::vector<double>> A3 = {
        {0.001, 1.0, 1.0},
        {1.0, 2.0, 1.0},
        {2.0, 1.0, 3.0}
    };
    std::vector<double> b3 = {2.001, 4.0, 6.0};
    std::vector<double> expected3 = {1.0, 1.0, 1.0};
    try {
        result = pivot_gauss(A3, b3);
        print_linear_result("Pivoted Gaussian", result, expected3, 1e-4);
    } catch (const std::exception& e) {
        std::cout << "  Pivoted Gaussian Error: " << e.what() << "\n";
    }

    // Test 4: Gauss-Jordan - When we need reduced row echelon form
    std::cout << "\n--- Test 4: Gauss-Jordan - 3x3 System (x=1, y=2, z=3) ---\n";
    std::vector<std::vector<double>> A4 = {
        {2.0, 1.0, 1.0},
        {1.0, 3.0, 2.0},
        {1.0, 1.0, 1.0}
    };
    std::vector<double> b4 = {7.0, 13.0, 6.0};  // Corrected: 1*1 + 3*2 + 2*3 = 13
    std::vector<double> expected4 = {1.0, 2.0, 3.0};
    try {
        result = gauss_jordan(A4, b4);
        print_linear_result("Gauss-Jordan", result, expected4);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Jordan Error: " << e.what() << "\n";
    }

    // Test 5: LU Decomposition - For systems solved multiple times
    std::cout << "\n--- Test 5: LU Decomposition - 4x4 System ---\n";
    std::vector<std::vector<double>> A5 = {
        {4.0, 3.0, 2.0, 1.0},
        {3.0, 4.0, 3.0, 2.0},
        {2.0, 3.0, 4.0, 3.0},
        {1.0, 2.0, 3.0, 4.0}
    };
    std::vector<double> b5 = {10.0, 12.0, 12.0, 10.0};
    std::vector<double> expected5 = {1.0, 1.0, 1.0, 1.0};
    try {
        result = lu_decomposition(A5, b5);
        print_linear_result("LU Decomposition", result, expected5);
    } catch (const std::exception& e) {
        std::cout << "  LU Decomposition Error: " << e.what() << "\n";
    }

    // Test 6: Cholesky - Symmetric Positive Definite matrix
    std::cout << "\n--- Test 6: Cholesky - SPD 4x4 System ---\n";
    std::vector<std::vector<double>> A6 = {
        {4.0, 2.0, 1.0, 0.5},
        {2.0, 5.0, 2.0, 1.0},
        {1.0, 2.0, 6.0, 2.0},
        {0.5, 1.0, 2.0, 7.0}
    };
    std::vector<double> b6 = {7.5, 10.0, 11.0, 10.5};
    std::vector<double> expected6 = {1.0, 1.0, 1.0, 1.0};
    try {
        result = cholesky_decomposition(A6, b6);
        print_linear_result("Cholesky Decomposition", result, expected6);
    } catch (const std::exception& e) {
        std::cout << "  Cholesky Decomposition Error: " << e.what() << "\n";
    }

    // Test 7: Thomas Algorithm - Large tridiagonal system
    std::cout << "\n--- Test 7: Thomas Algorithm - Tridiagonal 5x5 System ---\n";
    std::vector<std::vector<double>> A7 = {
        {2.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0, 0.0},
        {0.0, 1.0, 2.0, 1.0, 0.0},
        {0.0, 0.0, 1.0, 2.0, 1.0},
        {0.0, 0.0, 0.0, 1.0, 2.0}
    };
    std::vector<double> b7 = {3.0, 4.0, 4.0, 4.0, 3.0};
    std::vector<double> expected7 = {1.0, 1.0, 1.0, 1.0, 1.0};
    try {
        result = thomas_algorithm(A7, b7);
        print_linear_result("Thomas Algorithm", result, expected7);
    } catch (const std::exception& e) {
        std::cout << "  Thomas Algorithm Error: " << e.what() << "\n";
    }

    // Test 8: Jacobi - Diagonally dominant matrix
    std::cout << "\n--- Test 8: Jacobi - Diagonally Dominant 4x4 System ---\n";
    std::vector<std::vector<double>> A8 = {
        {10.0, 1.0, 2.0, 1.0},
        {1.0, 10.0, 1.0, 2.0},
        {2.0, 1.0, 10.0, 1.0},
        {1.0, 2.0, 1.0, 10.0}
    };
    std::vector<double> b8 = {14.0, 14.0, 14.0, 14.0};
    std::vector<double> expected8 = {1.0, 1.0, 1.0, 1.0};
    try {
        result = jacobi(A8, b8, max_iter, tol);
        print_linear_result("Jacobi", result, expected8);
    } catch (const std::exception& e) {
        std::cout << "  Jacobi Error: " << e.what() << "\n";
    }

    // Test 9: Gauss-Seidel - Similar to Jacobi but faster convergence
    std::cout << "\n--- Test 9: Gauss-Seidel - Diagonally Dominant 4x4 System ---\n";
    std::vector<std::vector<double>> A9 = {
        {10.0, 1.0, 1.0, 1.0},
        {1.0, 10.0, 1.0, 1.0},
        {1.0, 1.0, 10.0, 1.0},
        {1.0, 1.0, 1.0, 10.0}
    };
    std::vector<double> b9 = {13.0, 13.0, 13.0, 13.0};
    std::vector<double> expected9 = {1.0, 1.0, 1.0, 1.0};
    try {
        result = gauss_seidel(A9, b9, max_iter, tol);
        print_linear_result("Gauss-Seidel", result, expected9);
    } catch (const std::exception& e) {
        std::cout << "  Gauss-Seidel Error: " << e.what() << "\n";
    }

    // Test 10: SOR - System that needs acceleration
    std::cout << "\n--- Test 10: SOR - System Requiring Acceleration ---\n";
    std::vector<std::vector<double>> A10 = {
        {4.0, 1.0, 0.0},
        {1.0, 4.0, 1.0},
        {0.0, 1.0, 4.0}
    };
    std::vector<double> b10 = {5.0, 6.0, 5.0};
    std::vector<double> expected10 = {1.0, 1.0, 1.0};
    try {
        result = sor(A10, b10, max_iter, tol, 1.2);
        print_linear_result("SOR (omega=1.2)", result, expected10);
    } catch (const std::exception& e) {
        std::cout << "  SOR Error: " << e.what() << "\n";
    }

    // Test 11: Conjugate Gradient - Large SPD matrix
    std::cout << "\n--- Test 11: Conjugate Gradient - Large SPD 5x5 System ---\n";
    std::vector<std::vector<double>> A11 = {
        {5.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 5.0, 1.0, 0.0, 0.0},
        {0.0, 1.0, 5.0, 1.0, 0.0},
        {0.0, 0.0, 1.0, 5.0, 1.0},
        {0.0, 0.0, 0.0, 1.0, 5.0}
    };
    std::vector<double> b11 = {6.0, 7.0, 7.0, 7.0, 6.0};
    std::vector<double> expected11 = {1.0, 1.0, 1.0, 1.0, 1.0};
    try {
        result = conjugate_gradient(A11, b11, max_iter, tol);
        print_linear_result("Conjugate Gradient", result, expected11);
    } catch (const std::exception& e) {
        std::cout << "  Conjugate Gradient Error: " << e.what() << "\n";
    }
}

// ============================================
// INTERPOLATION TESTS
// ============================================

void print_interpolation_result(const std::string& method_name, const InterpolationResult& result, 
                                const std::vector<double>& expected, double tolerance = 1e-6) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Success: " << (result.success ? "Yes" : "No") << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    
    if (!result.success || result.interpolated_values.size() != expected.size()) {
        std::cout << "  ✗ FAILED (success=false or size mismatch)\n";
        return;
    }
    
    double max_error = 0.0;
    for (size_t i = 0; i < result.interpolated_values.size() && i < expected.size(); ++i) {
        double error = std::fabs(result.interpolated_values[i] - expected[i]);
        if (error > max_error) max_error = error;
    }
    std::cout << "  Max Error: " << max_error << "\n";
    if (result.success && max_error < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

void test_all_interpolation_methods() {
    std::cout << "\n\n========================================\n";
    std::cout << "TESTING ALL INTERPOLATION METHODS\n";
    std::cout << "========================================\n";

    // Test 1: Linear Interpolation - Random but valid data points (linear function)
    std::cout << "\n--- Test 1: Linear Interpolation - Linear Function f(x) = 2x + 1 ---\n";
    std::vector<double> x1 = {0.0, 1.5, 3.2, 4.8, 6.0};
    std::vector<double> y1;
    for (double xi : x1) y1.push_back(2.0 * xi + 1.0);
    std::vector<double> x_query1 = {0.75, 2.35, 4.0, 5.4};
    std::vector<double> expected1;
    for (double xq : x_query1) expected1.push_back(2.0 * xq + 1.0);
    
    try {
        InterpolationResult result = linear_interpolation(x1, y1, x_query1);
        print_interpolation_result("Linear Interpolation", result, expected1, 1e-6);
    } catch (const std::exception& e) {
        std::cout << "  Linear Interpolation Error: " << e.what() << "\n";
    }

    // Test 2: Lagrange Polynomial - Polynomial function f(x) = x^3
    std::cout << "\n--- Test 2: Lagrange Polynomial - Polynomial f(x) = x^3 ---\n";
    std::vector<double> x2 = {-2.5, -1.2, 0.3, 1.8, 3.1};
    std::vector<double> y2;
    for (double xi : x2) y2.push_back(xi * xi * xi);
    std::vector<double> x_query2 = {-1.85, 0.0, 1.05, 2.45};
    std::vector<double> expected2;
    for (double xq : x_query2) expected2.push_back(xq * xq * xq);
    
    try {
        InterpolationResult result = lagrange_polynomial(x2, y2, x_query2);
        print_interpolation_result("Lagrange Polynomial", result, expected2, 1e-3);
    } catch (const std::exception& e) {
        std::cout << "  Lagrange Polynomial Error: " << e.what() << "\n";
    }

    // Test 3: Newton Divided Difference - Exponential function
    std::cout << "\n--- Test 3: Newton Divided Difference - Exponential f(x) = exp(0.5x) ---\n";
    std::vector<double> x3 = {0.2, 0.8, 1.5, 2.3, 3.1};
    std::vector<double> y3;
    for (double xi : x3) y3.push_back(std::exp(0.5 * xi));
    std::vector<double> x_query3 = {0.5, 1.15, 1.9, 2.7};
    std::vector<double> expected3;
    for (double xq : x_query3) expected3.push_back(std::exp(0.5 * xq));
    
    try {
        InterpolationResult result = newton_divided_difference(x3, y3, x_query3);
        print_interpolation_result("Newton Divided Difference", result, expected3, 1e-3);
    } catch (const std::exception& e) {
        std::cout << "  Newton Divided Difference Error: " << e.what() << "\n";
    }

    // Test 4: Cubic Spline - Smooth function f(x) = sin(x) + 0.5x
    std::cout << "\n--- Test 4: Cubic Spline - Smooth Function f(x) = sin(x) + 0.5x ---\n";
    std::vector<double> x4 = {0.0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2};
    std::vector<double> y4;
    for (double xi : x4) y4.push_back(std::sin(xi) + 0.5 * xi);
    std::vector<double> x_query4 = {0.35, 1.05, 1.75, 2.45, 3.15};
    std::vector<double> expected4;
    for (double xq : x_query4) expected4.push_back(std::sin(xq) + 0.5 * xq);
    
    try {
        InterpolationResult result = cubic_spline(x4, y4, x_query4);
        print_interpolation_result("Cubic Spline (Natural)", result, expected4, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  Cubic Spline Error: " << e.what() << "\n";
    }
}

// ============================================
// ODE SOLVER TESTS
// ============================================

void print_ode_result(const std::string& method_name, const ODEResult& result, 
                      double expected_final, double tolerance = 1e-3) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Steps: " << result.steps << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    std::cout << "  Max Error: " << result.max_error << "\n";
    if (!result.solution_single.empty()) {
        double final_value = result.solution_single.back();
        std::cout << "  Final Value: " << final_value << " (expected: " << expected_final << ")\n";
        double error = std::fabs(final_value - expected_final);
        std::cout << "  Error: " << error << "\n";
        if (result.converged && error < tolerance) {
            std::cout << "  ✓ PASSED\n";
        } else {
            std::cout << "  ✗ FAILED\n";
        }
    } else {
        std::cout << "  ✗ FAILED (no solution)\n";
    }
}

void test_all_ode_solver_methods() {
    std::cout << "\n\n========================================\n";
    std::cout << "TESTING ALL ODE SOLVER METHODS\n";
    std::cout << "========================================\n";

    // Test 1: Euler - Simple exponential decay: dy/dt = -y, y(0) = 1, solution: y(t) = exp(-t)
    std::cout << "\n--- Test 1: Euler - Exponential Decay (dy/dt = -y, y(0)=1) ---\n";
    ODEFunc f1 = [](double t, double y) { (void)t; return -y; };
    double t0 = 0.0, y0 = 1.0, t_end = 1.0, h = 0.01;
    double expected1 = std::exp(-t_end);
    
    try {
        ODEResult result = euler(f1, t0, y0, t_end, h);
        print_ode_result("Euler", result, expected1, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  Euler Error: " << e.what() << "\n";
    }

    // Test 2: RK4 - Same problem for comparison
    std::cout << "\n--- Test 2: RK4 - Exponential Decay (dy/dt = -y, y(0)=1) ---\n";
    try {
        ODEResult result = rk4(f1, t0, y0, t_end, h);
        print_ode_result("RK4", result, expected1, 1e-4);
    } catch (const std::exception& e) {
        std::cout << "  RK4 Error: " << e.what() << "\n";
    }

    // Test 3: RK45 - Adaptive method
    std::cout << "\n--- Test 3: RK45 - Adaptive Exponential Decay ---\n";
    try {
        ODEResult result = rk45(f1, t0, y0, t_end, 0.1, 1e-6);
        print_ode_result("RK45", result, expected1, 1e-4);
    } catch (const std::exception& e) {
        std::cout << "  RK45 Error: " << e.what() << "\n";
    }

    // Test 4: RK78 - High-order adaptive method
    std::cout << "\n--- Test 4: RK78 - High-Order Adaptive Method ---\n";
    try {
        ODEResult result = rk78(f1, t0, y0, t_end, 0.1, 1e-6);
        print_ode_result("RK78", result, expected1, 1e-4);
    } catch (const std::exception& e) {
        std::cout << "  RK78 Error: " << e.what() << "\n";
    }

    // Test 5: BDF - Stiff ODE: dy/dt = -100y, y(0) = 1, solution: y(t) = exp(-100t)
    std::cout << "\n--- Test 5: BDF - Stiff ODE (dy/dt = -100y, y(0)=1) ---\n";
    ODEFunc f2 = [](double t, double y) { (void)t; return -100.0 * y; };
    double expected2 = std::exp(-100.0 * 0.1);
    try {
        ODEResult result = bdf(f2, t0, y0, 0.1, 0.001, 2);
        print_ode_result("BDF", result, expected2, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  BDF Error: " << e.what() << "\n";
    }

    // Test 6: Radau - Implicit method for stiff systems
    std::cout << "\n--- Test 6: Radau - Implicit Method for Stiff ODE ---\n";
    try {
        ODEResult result = radau(f2, t0, y0, 0.1, 0.001, 1e-6);
        print_ode_result("Radau", result, expected2, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  Radau Error: " << e.what() << "\n";
    }
}

// ============================================
// OPTIMIZATION TESTS
// ============================================

void print_optimization_result(const std::string& method_name, const OptimizationResult& result, 
                               const std::vector<double>& expected, double tolerance = 1e-3) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << method_name << ":\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  CPU Time: " << std::scientific << result.cpu_time_sec << " sec\n";
    std::cout << "  FLOP Count: " << result.flop_count << "\n";
    std::cout << "  Optimal Value: " << result.optimal_value << "\n";
    
    double max_error = 0.0;
    for (size_t i = 0; i < result.optimal_point.size() && i < expected.size(); ++i) {
        double error = std::fabs(result.optimal_point[i] - expected[i]);
        if (error > max_error) max_error = error;
    }
    std::cout << "  Max Error: " << max_error << "\n";
    if (result.converged && max_error < tolerance) {
        std::cout << "  ✓ PASSED\n";
    } else {
        std::cout << "  ✗ FAILED\n";
    }
}

void test_all_optimization_methods() {
    std::cout << "\n\n========================================\n";
    std::cout << "TESTING ALL OPTIMIZATION METHODS\n";
    std::cout << "========================================\n";

    // Test 1: Steepest Descent - Simple quadratic: f(x) = sum(x_i^2)
    std::cout << "\n--- Test 1: Steepest Descent - f(x) = sum(x_i^2) ---\n";
    ObjectiveFunc f1 = [](const std::vector<double>& x) {
        double sum = 0.0;
        for (double xi : x) sum += xi * xi;
        return sum;
    };
    GradientFunc grad_f1 = [](const std::vector<double>& x) {
        std::vector<double> grad(x.size());
        for (size_t i = 0; i < x.size(); ++i) grad[i] = 2.0 * x[i];
        return grad;
    };
    std::vector<double> x0_1 = {5.0, 5.0, 5.0};
    std::vector<double> expected1 = {0.0, 0.0, 0.0};
    
    try {
        OptimizationResult result = steepest_descent(f1, grad_f1, x0_1, 0.1, 1e-6, 1000);
        print_optimization_result("Steepest Descent", result, expected1, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  Steepest Descent Error: " << e.what() << "\n";
    }

    // Test 2: BFGS - Simple quadratic (easier than Rosenbrock)
    std::cout << "\n--- Test 2: BFGS - Simple Quadratic f(x) = sum(x_i^2) ---\n";
    std::vector<double> x0_2 = {3.0, 3.0};
    std::vector<double> expected2 = {0.0, 0.0};
    
    try {
        OptimizationResult result = bfgs(f1, grad_f1, x0_2, 1e-6, 1000);
        print_optimization_result("BFGS", result, expected2, 1e-2);
    } catch (const std::exception& e) {
        std::cout << "  BFGS Error: " << e.what() << "\n";
    }

    // Test 3: Genetic Algorithm - Multi-modal function
    std::cout << "\n--- Test 3: Genetic Algorithm - f(x) = sum(x_i^2) ---\n";
    std::vector<double> lower_bounds = {-5.0, -5.0, -5.0};
    std::vector<double> upper_bounds = {5.0, 5.0, 5.0};
    
    try {
        OptimizationResult result = genetic_algorithm(f1, 3, lower_bounds, upper_bounds, 50, 100, 0.1, 0.8);
        print_optimization_result("Genetic Algorithm", result, expected1, 1e-1);
    } catch (const std::exception& e) {
        std::cout << "  Genetic Algorithm Error: " << e.what() << "\n";
    }

    // Test 4: PSO - Particle Swarm Optimization
    std::cout << "\n--- Test 4: PSO - f(x) = sum(x_i^2) ---\n";
    try {
        OptimizationResult result = pso(f1, 3, lower_bounds, upper_bounds, 30, 100, 0.7, 1.5, 1.5);
        print_optimization_result("PSO", result, expected1, 1e-1);
    } catch (const std::exception& e) {
        std::cout << "  PSO Error: " << e.what() << "\n";
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

    // Test all interpolation methods
    test_all_interpolation_methods();

    // Test all ODE solver methods
    test_all_ode_solver_methods();

    // Test all optimization methods
    test_all_optimization_methods();

    std::cout << "\n\n========================================\n";
    std::cout << "TESTING COMPLETE\n";
    std::cout << "========================================\n";

    return 0;
}
