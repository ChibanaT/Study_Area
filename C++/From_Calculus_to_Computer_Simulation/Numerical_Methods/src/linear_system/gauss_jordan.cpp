#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace numerical {

//----------------------------
// Gauss-Jordan Elimination
//----------------------------
LinearResult gauss_jordan(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    // Create augmented matrix
    std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) aug[i][j] = A_in[i][j];
        aug[i][n] = b_in[i];
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        int maxRow = i;
        double maxVal = std::fabs(aug[i][i]);
        for (int k = i+1; k < n; ++k) {
            if (std::fabs(aug[k][i]) > maxVal) {
                maxVal = std::fabs(aug[k][i]);
                maxRow = k;
            }
        }
        if (maxRow != i) std::swap(aug[i], aug[maxRow]);

        if (std::fabs(aug[i][i]) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in Gauss-Jordan");

        // Normalize pivot row
        double pivot = aug[i][i];
        for (int j = i; j <= n; ++j) {
            aug[i][j] /= pivot;
            result.flop_count += 1; // division
        }

        // Eliminate other rows
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = aug[k][i];
            for (int j = i; j <= n; ++j) {
                aug[k][j] -= factor * aug[i][j];
                result.flop_count += 2; // multiply + subtract
            }
        }
    }

    // Extract solution
    result.solution.resize(n);
    for (int i = 0; i < n; ++i) result.solution[i] = aug[i][n];

    result.iterations = 1;  // Direct method
    result.converged = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
