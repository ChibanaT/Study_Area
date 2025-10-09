#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult newton_raphson(const Func &f, const Func &df, double x0, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double x = x0;
    double fx = f(x);
    double dfx = df(x);
    result.flop_count += 2;

    for (int k = 0; k < max_iter; ++k) {
        if (std::fabs(dfx) < 1e-14)
            throw std::runtime_error("Derivative too small — division by zero risk");

        // Newton-Raphson update
        double x_next = x - fx / dfx;
        double fx_next = f(x_next);
        result.history.push_back(x_next);
        ++result.iterations;
        result.flop_count += 5; // division, subtraction, and new function eval

        // Check convergence
        if (std::fabs(fx_next) < tol || std::fabs(x_next - x) < tol) {
            result.root = x_next;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Update for next iteration
        x = x_next;
        fx = fx_next;
        dfx = df(x);
        result.flop_count += 1; // derivative eval
    }

    // If not converged
    result.root = x;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
