#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace numerical {

//----------------------------
// Gaussian Elimination with Partial Pivoting
//----------------------------
LinearResult pivot_gauss(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    // Copy inputs to work on
    std::vector<std::vector<double>> A = A_in;
    std::vector<double> b = b_in;
    result.solution.resize(n);

    // Forward elimination with partial pivoting
    for (int k = 0; k < n-1; ++k) {
        // Find the row with the largest pivot
        int maxRow = k;
        double maxVal = std::fabs(A[k][k]);
        for (int i = k+1; i < n; ++i) {
            if (std::fabs(A[i][k]) > maxVal) {
                maxVal = std::fabs(A[i][k]);
                maxRow = i;
            }
        }

        // Swap rows if necessary
        if (maxRow != k) {
            std::swap(A[k], A[maxRow]);
            std::swap(b[k], b[maxRow]);
        }

        if (std::fabs(A[k][k]) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in pivoted Gaussian elimination");

        // Eliminate below pivot
        for (int i = k+1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            result.flop_count += 1; // division
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
                result.flop_count += 2; // multiply + subtract
            }
            b[i] -= factor * b[k];
            result.flop_count += 2; // multiply + subtract
        }
    }

    // Back substitution
    for (int i = n-1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j < n; ++j) {
            sum += A[i][j] * result.solution[j];
            result.flop_count += 2; // multiply + add
        }
        result.solution[i] = (b[i] - sum) / A[i][i];
        result.flop_count += 1; // division
    }

    result.iterations = 1; // Direct method
    result.converged = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
