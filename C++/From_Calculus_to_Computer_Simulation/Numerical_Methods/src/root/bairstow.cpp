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
        throw std::invalid_argument("Polynomial degree must be >= 2 for Bairstow");

    std::vector<double> a = coeffs; // working copy
    std::vector<double> b(n+1);
    double dr = 0.0, ds = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Synthetic division to compute b coefficients
        b[n] = a[n];
        b[n-1] = a[n-1] + r*b[n];
        for (int i = n-2; i >= 0; --i) {
            b[i] = a[i] + r*b[i+1] + s*b[i+2];
        }
        result.flop_count += 4 * n; // rough estimate of operations

        // Compute corrections dr, ds
        double det = b[2]*b[2] - b[3]*b[1];
        if (std::fabs(det) < 1e-14)
            throw std::runtime_error("Determinant too small in Bairstow iteration");

        dr = (-b[1]*b[2] + b[0]*b[3]) / det;
        ds = (-b[0]*b[2] + b[1]*b[1]) / det;

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
