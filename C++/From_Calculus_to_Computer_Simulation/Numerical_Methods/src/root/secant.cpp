#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult secant(const Func &f, double a, double b, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double fa = f(a);
    double fb = f(b);
    result.flop_count += 2;

    if (std::fabs(fa) < tol) {
        result.root = a;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (std::fabs(fb) < tol) {
        result.root = b;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }

    double c = b;
    double fc = fb;

    for (int k = 0; k < max_iter; ++k) {
        if (std::fabs(fb - fa) < 1e-14)
            throw std::runtime_error("Division by zero in secant method");

        // Secant formula
        c = b - fb * (b - a) / (fb - fa);
        fc = f(c);
        result.history.push_back(c);
        ++result.iterations;
        result.flop_count += 6; // multiplications, subtractions, division, and function call

        // Check convergence
        if (std::fabs(fc) < tol || std::fabs(c - b) < tol) {
            result.root = c;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Shift points for next iteration
        a = b;
        fa = fb;
        b = c;
        fb = fc;
    }

    // If not converged
    result.root = c;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
