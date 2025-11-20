#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

//----------------------------
// Jacobi Iterative Method
//----------------------------
LinearResult jacobi(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in, int max_iter, double tol) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    // Initialize solution vector
    result.solution.resize(n, 0.0);
    std::vector<double> x_old(n, 0.0);
    std::vector<double> x_new(n);

    // Check diagonal dominance (warning only, not required)
    for (int i = 0; i < n; ++i) {
        if (std::fabs(A_in[i][i]) < 1e-12)
            throw std::runtime_error("Zero diagonal element encountered in Jacobi method");
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute new approximation
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A_in[i][j] * x_old[j];
                    result.flop_count += 2; // multiply + add
                }
            }
            x_new[i] = (b_in[i] - sum) / A_in[i][i];
            result.flop_count += 2; // subtraction + division
        }

        ++result.iterations;

        // Check convergence: ||x_new - x_old|| < tol
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = std::fabs(x_new[i] - x_old[i]);
            if (diff > error) error = diff;
            result.flop_count += 1; // subtraction and comparison
        }

        if (error < tol) {
            result.solution = x_new;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Update for next iteration
        x_old = x_new;
    }

    // If not converged
    result.solution = x_new;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical

