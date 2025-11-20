#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

namespace numerical {

// Helper function: evaluates polynomial using Horner's method
inline double horner(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    for (auto c : coeffs) result = result * x + c;
    return result;
}

RootResult bairstow(const std::vector<double>& coeffs, double r, double s, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = coeffs.size() - 1;
    if (n < 2)
        throw std::invalid_argument("Polynomial degree must be >= 2 for Bairstow method");

    std::vector<double> a = coeffs; // working copy
    std::vector<double> b(n+1);
    std::vector<double> c(n+1);
    double dr = 0.0, ds = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Synthetic division to compute b coefficients
        b[n] = a[n];
        b[n-1] = a[n-1] + r*b[n];
        for (int i = n-2; i >= 0; --i) {
            b[i] = a[i] + r*b[i+1] + s*b[i+2];
        }
        result.flop_count += 4 * n; // rough estimate of operations

        // Compute c coefficients (derivatives with respect to r and s)
        c[n] = 0.0;
        c[n-1] = b[n];
        for (int i = n-2; i >= 0; --i) {
            c[i] = b[i+1] + r*c[i+1] + s*c[i+2];
        }
        result.flop_count += 3 * n; // rough estimate

        // Compute corrections dr, ds using Newton's method
        // Solve: c[1]*dr + c[2]*ds = -b[0]
        //        c[2]*dr + c[3]*ds = -b[1]
        // For n=2, we need to handle the case where c[3] might not exist
        // Use c[2] as fallback or ensure proper indexing
        double c1 = (n >= 1) ? c[1] : 0.0;
        double c2 = (n >= 2) ? c[2] : 0.0;
        double c3 = (n >= 3) ? c[3] : c[2]; // For n=2, use c[2] as approximation
        
        double det = c1*c3 - c2*c2;
        if (std::fabs(det) < 1e-14)
            throw std::runtime_error("Determinant too small in Bairstow iteration");

        dr = (-b[0]*c3 + b[1]*c2) / det;
        ds = (-b[1]*c1 + b[0]*c2) / det;
        result.flop_count += 8; // multiplications, additions, divisions

        r += dr;
        s += ds;
        ++result.iterations;

        // Check convergence
        if (std::fabs(dr) < tol && std::fabs(ds) < tol) {
            // Roots of quadratic factor x^2 - r x - s = 0
            double disc = r*r + 4*s;
            if (disc >= 0) {
                result.root = 0.5 * (r + std::sqrt(disc));
                result.history.push_back(result.root);
                result.history.push_back(0.5 * (r - std::sqrt(disc)));
            } else {
                // complex roots
                result.root = r / 2.0; // return real part of first root
                result.history.push_back(r / 2.0 + std::sqrt(-disc)/2.0); // complex part in history
                result.history.push_back(r / 2.0 - std::sqrt(-disc)/2.0);
            }

            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            result.cpu_time_sec = std::chrono::duration<double>(end_time - start_time).count();
            return result;
        }
    }

    // If not converged
    result.root = r; // approximate
    auto end_time = std::chrono::high_resolution_clock::now();
    result.cpu_time_sec = std::chrono::duration<double>(end_time - start_time).count();
    return result;
}

} // namespace numerical
