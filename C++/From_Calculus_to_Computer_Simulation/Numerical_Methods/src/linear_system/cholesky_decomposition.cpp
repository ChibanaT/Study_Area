#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

//----------------------------
// Cholesky Decomposition (for symmetric positive definite matrices)
//----------------------------
LinearResult cholesky_decomposition(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    // Check symmetry (approximate)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (std::fabs(A_in[i][j] - A_in[j][i]) > 1e-10)
                throw std::invalid_argument("Matrix A must be symmetric for Cholesky decomposition");
        }
    }

    // Cholesky decomposition: A = L * L^T
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<double> y(n);

    // Compute L
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k < i; ++k) {
            sum += L[i][k] * L[i][k];
            result.flop_count += 2; // multiply + add
        }
        double diag = A_in[i][i] - sum;
        if (diag <= 0.0)
            throw std::runtime_error("Matrix is not positive definite");
        L[i][i] = std::sqrt(diag);
        result.flop_count += 1; // sqrt

        for (int j = i + 1; j < n; ++j) {
            sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += L[j][k] * L[i][k];
                result.flop_count += 2; // multiply + add
            }
            L[j][i] = (A_in[j][i] - sum) / L[i][i];
            result.flop_count += 2; // subtraction + division
        }
    }

    // Forward substitution: Ly = b
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
            result.flop_count += 2; // multiply + add
        }
        y[i] = (b_in[i] - sum) / L[i][i];
        result.flop_count += 2; // subtraction + division
    }

    // Back substitution: L^T x = y
    result.solution.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += L[j][i] * result.solution[j];
            result.flop_count += 2; // multiply + add
        }
        result.solution[i] = (y[i] - sum) / L[i][i];
        result.flop_count += 2; // subtraction + division
    }

    result.iterations = 1;  // Direct method
    result.converged = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical

