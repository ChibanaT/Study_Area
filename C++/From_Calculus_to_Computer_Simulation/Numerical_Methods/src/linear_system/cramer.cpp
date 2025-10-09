#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>

namespace numerical {

// Compute determinant recursively 
double determinant(const std::vector<std::vector<double>>& mat) {
    int n = mat.size();
    if (n == 1) return mat[0][0];
    if (n == 2) return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];

    double det = 0.0;
    for (int p = 0; p < n; ++p) {
        std::vector<std::vector<double>> submat(n-1, std::vector<double>(n-1));
        for (int i = 1; i < n; ++i) {
            int colIndex = 0;
            for (int j = 0; j < n; ++j) {
                if (j == p) continue;
                submat[i-1][colIndex] = mat[i][j];
                ++colIndex;
            }
        }
        det += (p % 2 == 0 ? 1 : -1) * mat[0][p] * determinant(submat);
    }
    return det;
}

//----------------------------
// Cramer's Rule
//----------------------------
LinearResult cramer(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 1;      // Single computation
    result.converged = false;
    result.flop_count = 0;

    int n = A.size();
    if (n == 0 || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    double detA = determinant(A);
    result.flop_count += n * n; // Rough estimate
    if (detA == 0.0) {
        result.solution = std::vector<double>(n, 0.0);
        result.converged = false;
    } else {
        result.solution.resize(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::vector<double>> Ai = A;
            for (int j = 0; j < n; ++j) Ai[j][i] = b[j];
            result.solution[i] = determinant(Ai) / detA;
            result.flop_count += n * n; // Rough estimate
        }
        result.converged = true;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
