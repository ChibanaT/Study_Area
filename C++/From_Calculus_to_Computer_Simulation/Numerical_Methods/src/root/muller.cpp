#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult muller(const Func &f, double x0, double x1, double x2, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double fx0 = f(x0);
    double fx1 = f(x1);
    double fx2 = f(x2);
    result.flop_count += 3;

    // Check if any initial point is a root
    if (std::fabs(fx0) < tol) {
        result.root = x0;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (std::fabs(fx1) < tol) {
        result.root = x1;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (std::fabs(fx2) < tol) {
        result.root = x2;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }

    double x = x2;
    double fx = fx2;

    for (int k = 0; k < max_iter; ++k) {
        // Compute divided differences
        double h0 = x1 - x0;
        double h1 = x2 - x1;
        
        if (std::fabs(h0) < 1e-14 || std::fabs(h1) < 1e-14)
            throw std::runtime_error("Division by zero in Muller method: points too close");
        
        double d0 = (fx1 - fx0) / h0;
        double d1 = (fx2 - fx1) / h1;
        double a = (d1 - d0) / (h1 + h0);
        double b = a * h1 + d1;
        double c = fx2;
        result.flop_count += 10; // divisions, subtractions, multiplications

        // Compute discriminant
        double discriminant = b * b - 4.0 * a * c;
        result.flop_count += 3; // multiplications and subtraction

        if (std::fabs(a) < 1e-14) {
            // Linear case: use secant method
            if (std::fabs(fx2 - fx1) < 1e-14)
                throw std::runtime_error("Division by zero in Muller method");
            x = x2 - fx2 * (x2 - x1) / (fx2 - fx1);
            result.flop_count += 4;
        } else if (discriminant < 0) {
            // Complex roots - use real part
            x = x2 - 2.0 * c / (b + std::sqrt(-discriminant) * (b >= 0 ? 1 : -1));
            result.flop_count += 4;
        } else {
            // Two real roots - choose the one closer to x2
            double sqrt_disc = std::sqrt(discriminant);
            double denom1 = b + sqrt_disc;
            double denom2 = b - sqrt_disc;
            double denom = (std::fabs(denom1) > std::fabs(denom2)) ? denom1 : denom2;
            x = x2 - 2.0 * c / denom;
            result.flop_count += 5;
        }

        fx = f(x);
        result.history.push_back(x);
        ++result.iterations;
        result.flop_count += 1; // function evaluation

        // Check convergence
        if (std::fabs(fx) < tol || std::fabs(x - x2) < tol) {
            result.root = x;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Shift points for next iteration
        x0 = x1;
        fx0 = fx1;
        x1 = x2;
        fx1 = fx2;
        x2 = x;
        fx2 = fx;
    }

    // If not converged
    result.root = x;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical

