#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult modified_false_position(const Func &f, double a, double b, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double fa = f(a), fb = f(b);
    result.flop_count += 2;

    if (fa == 0.0) {
        result.root = a;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (fb == 0.0) {
        result.root = b;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (fa * fb > 0.0)
        throw std::invalid_argument("f(a) and f(b) must have opposite signs");

    double c = a;
    double fc = fa;
    bool modified_a = false;
    bool modified_b = false;

    for (int k = 0; k < max_iter; ++k) {
        // Regula Falsi formula
        c = (a * fb - b * fa) / (fb - fa);
        fc = f(c);
        result.history.push_back(c);
        ++result.iterations;
        result.flop_count += 6; // multiplications/divisions and one function call

        // Check convergence
        if (std::fabs(fc) < tol || std::fabs(b - a) < tol) {
            result.root = c;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Determine interval update
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
            if (modified_a) fa *= 0.5; // scale stagnant endpoint
            modified_a = false;
            modified_b = true;
        } else {
            a = c;
            fa = fc;
            if (modified_b) fb *= 0.5;
            modified_b = false;
            modified_a = true;
        }

        result.flop_count += 2; // includes scaling operation and comparison
    }

    // If not converged
    result.root = c;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
