#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult halley(const Func &f, const Func &df, const Func &d2f, double x0, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double x = x0;
    double fx = f(x);
    double dfx = df(x);
    double ddfx = d2f(x);
    result.flop_count += 3;

    for (int k = 0; k < max_iter; ++k) {
        if (std::fabs(dfx) < 1e-14)
            throw std::runtime_error("Derivative too small — division by zero risk");

        // Halley update:
        // x_{n+1} = x_n - (2 f f') / (2 (f')^2 - f f'')
        double numerator = 2.0 * fx * dfx;
        double denominator = 2.0 * dfx * dfx - fx * ddfx;

        if (std::fabs(denominator) < 1e-14)
            throw std::runtime_error("Denominator too small — numerical instability");

        double x_next = x - numerator / denominator;
        double fx_next = f(x_next);
        result.history.push_back(x_next);
        ++result.iterations;
        result.flop_count += 8; // function evaluations + arithmetic ops

        // Convergence check
        if (std::fabs(fx_next) < tol || std::fabs(x_next - x) < tol) {
            result.root = x_next;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Prepare for next iteration
        x = x_next;
        fx = fx_next;
        dfx = df(x);
        ddfx = d2f(x);
        result.flop_count += 2; // derivative evaluations
    }

    // If not converged
    result.root = x;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
