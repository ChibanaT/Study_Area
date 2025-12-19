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
    std::vector<double> b(n+1, 0.0);
    std::vector<double> c(n+1, 0.0);
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
        // The system to solve is:
        // c[n-1]*dr + c[n-2]*ds = -b[0]
        // c[n-2]*dr + c[n-3]*ds = -b[1]
        // But we need to use the correct indices for the derivatives
        double c1 = (n >= 1) ? c[n-1] : 0.0;  // dc/dr
        double c2 = (n >= 2) ? c[n-2] : 0.0;  // dc/ds  
        double c3 = (n >= 2) ? c[n-2] : 0.0;  // same as c2 for second equation
        double c4 = (n >= 3) ? c[n-3] : 0.0;  // second derivative term
        
        // Solve the 2x2 system: [c1 c2] [dr]   [-b[0]]
        //                       [c2 c4] [ds] = [-b[1]]
        double det = c1*c4 - c2*c2;
        if (std::fabs(det) < 1e-14) {
            // If determinant is too small, use a simpler update
            dr = -b[0] * 0.1;
            ds = -b[1] * 0.1;
        } else {
            dr = (-b[0]*c4 + b[1]*c2) / det;
            ds = (-b[1]*c1 + b[0]*c2) / det;
        }
        result.flop_count += 8; // multiplications, additions, divisions

        r += dr;
        s += ds;
        ++result.iterations;

        // Check convergence
        if (std::fabs(dr) < tol && std::fabs(ds) < tol) {
            // Roots of quadratic factor x^2 - r x - s = 0
            // Using quadratic formula: x = (r ± sqrt(r^2 + 4s)) / 2
            double disc = r*r + 4*s;
            if (disc >= 0) {
                double sqrt_disc = std::sqrt(disc);
                double root1 = 0.5 * (r + sqrt_disc);
                double root2 = 0.5 * (r - sqrt_disc);
                result.root = root1;
                result.history.push_back(root1);
                result.history.push_back(root2);
            } else {
                // complex roots - return real part
                result.root = r / 2.0;
                result.history.push_back(r / 2.0);
                result.history.push_back(r / 2.0);
            }

            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            result.cpu_time_sec = std::chrono::duration<double>(end_time - start_time).count();
            return result;
        }
        
        // Prevent divergence
        if (std::fabs(dr) > 10.0 || std::fabs(ds) > 10.0 || std::isnan(dr) || std::isnan(ds)) {
            break;
        }
    }

    // If not converged, try to extract a root from the last iteration
    // Roots of quadratic factor x^2 - r x - s = 0
    double disc = r*r + 4*s;
    if (disc >= 0 && !std::isnan(r) && !std::isnan(s)) {
        double sqrt_disc = std::sqrt(disc);
        result.root = 0.5 * (r + sqrt_disc);
        result.history.push_back(result.root);
        result.history.push_back(0.5 * (r - sqrt_disc));
        result.converged = false; // Mark as not converged but return best guess
    } else {
        result.root = 0.0;
        result.converged = false;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.cpu_time_sec = std::chrono::duration<double>(end_time - start_time).count();
    return result;
}

} // namespace numerical
