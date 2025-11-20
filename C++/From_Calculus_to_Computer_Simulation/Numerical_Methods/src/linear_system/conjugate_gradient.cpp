#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

//----------------------------
// Conjugate Gradient Method (for symmetric positive definite matrices)
//----------------------------
LinearResult conjugate_gradient(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in, int max_iter, double tol) {
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
    std::vector<double> r(n);  // residual
    std::vector<double> p(n);  // search direction
    std::vector<double> Ap(n); // A * p

    // Compute initial residual: r = b - A * x
    for (int i = 0; i < n; ++i) {
        double Ax = 0.0;
        for (int j = 0; j < n; ++j) {
            Ax += A_in[i][j] * result.solution[j];
            result.flop_count += 2; // multiply + add
        }
        r[i] = b_in[i] - Ax;
        result.flop_count += 1; // subtraction
    }

    // Initial search direction
    p = r;
    double r_dot_r = 0.0;
    for (int i = 0; i < n; ++i) {
        r_dot_r += r[i] * r[i];
        result.flop_count += 2; // multiply + add
    }

    // Check initial convergence
    if (std::sqrt(r_dot_r) < tol) {
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute Ap = A * p
        for (int i = 0; i < n; ++i) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; ++j) {
                Ap[i] += A_in[i][j] * p[j];
                result.flop_count += 2; // multiply + add
            }
        }

        // Compute p^T * Ap
        double pAp = 0.0;
        for (int i = 0; i < n; ++i) {
            pAp += p[i] * Ap[i];
            result.flop_count += 2; // multiply + add
        }

        if (std::fabs(pAp) < 1e-14)
            throw std::runtime_error("Division by zero in conjugate gradient");

        // Compute step size: alpha = (r^T * r) / (p^T * Ap)
        double alpha = r_dot_r / pAp;
        result.flop_count += 1; // division

        // Update solution: x = x + alpha * p
        for (int i = 0; i < n; ++i) {
            result.solution[i] += alpha * p[i];
            result.flop_count += 2; // multiply + add
        }

        // Update residual: r = r - alpha * Ap
        double r_dot_r_new = 0.0;
        for (int i = 0; i < n; ++i) {
            r[i] -= alpha * Ap[i];
            r_dot_r_new += r[i] * r[i];
            result.flop_count += 3; // multiply, subtract, multiply, add
        }

        ++result.iterations;

        // Check convergence
        if (std::sqrt(r_dot_r_new) < tol) {
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Compute beta for next search direction: beta = (r_new^T * r_new) / (r^T * r)
        double beta = r_dot_r_new / r_dot_r;
        result.flop_count += 1; // division

        // Update search direction: p = r + beta * p
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
            result.flop_count += 2; // multiply + add
        }

        r_dot_r = r_dot_r_new;
    }

    // If not converged
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical

