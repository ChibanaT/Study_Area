#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

//----------------------------
// LU Decomposition (Doolittle's method)
//----------------------------
LinearResult lu_decomposition(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<double> y(n);

    // Copy A
    std::vector<std::vector<double>> A = A_in;

    // Decomposition
    for (int i = 0; i < n; ++i) {
        // Upper triangular
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += L[i][j] * U[j][k];
                result.flop_count += 2; // multiply + add
            }
            U[i][k] = A[i][k] - sum;
            result.flop_count += 1; // subtraction
        }

        if (std::fabs(U[i][i]) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in LU decomposition");

        // Lower triangular
        L[i][i] = 1.0;
        for (int k = i+1; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += L[k][j] * U[j][i];
                result.flop_count += 2; // multiply + add
            }
            L[k][i] = (A[k][i] - sum) / U[i][i];
            result.flop_count += 1; // division
        }
    }

    // Forward substitution Ly = b
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
            result.flop_count += 2; // multiply + add
        }
        y[i] = b_in[i] - sum;
        result.flop_count += 1; // subtraction
    }

    // Back substitution Ux = y
    result.solution.resize(n);
    for (int i = n-1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j < n; ++j) {
            sum += U[i][j] * result.solution[j];
            result.flop_count += 2; // multiply + add
        }
        result.solution[i] = (y[i] - sum) / U[i][i];
        result.flop_count += 1; // division
    }

    result.iterations = 1;  // Direct method
    result.converged = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
