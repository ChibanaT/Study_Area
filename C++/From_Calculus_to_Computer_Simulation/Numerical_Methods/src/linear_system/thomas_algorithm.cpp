#include "../../include/linear_system.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

//----------------------------
// Thomas Algorithm (for tridiagonal systems)
//----------------------------
LinearResult thomas_algorithm(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in) {
    auto start_time = std::chrono::high_resolution_clock::now();

    LinearResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    int n = A_in.size();
    if (n == 0 || A_in[0].size() != n || b_in.size() != n)
        throw std::invalid_argument("Matrix A must be square and vector b compatible");

    // Extract tridiagonal elements
    // For a tridiagonal matrix, we assume:
    // A[i][i-1] = lower diagonal (a_i), i = 1..n-1
    // A[i][i]   = main diagonal (b_i), i = 0..n-1
    // A[i][i+1] = upper diagonal (c_i), i = 0..n-2
    std::vector<double> a(n, 0.0);  // lower diagonal
    std::vector<double> b(n);       // main diagonal
    std::vector<double> c(n, 0.0);  // upper diagonal
    std::vector<double> d = b_in;   // right-hand side

    // Extract diagonal elements
    for (int i = 0; i < n; ++i) {
        b[i] = A_in[i][i];
        if (i > 0) a[i] = A_in[i][i-1];
        if (i < n-1) c[i] = A_in[i][i+1];
    }

    // Forward elimination
    for (int i = 1; i < n; ++i) {
        if (std::fabs(b[i-1]) < 1e-12)
            throw std::runtime_error("Zero pivot encountered in Thomas algorithm");

        double factor = a[i] / b[i-1];
        b[i] -= factor * c[i-1];
        d[i] -= factor * d[i-1];
        result.flop_count += 5; // division, 2 multiplications, 2 subtractions
    }

    // Back substitution
    result.solution.resize(n);
    if (std::fabs(b[n-1]) < 1e-12)
        throw std::runtime_error("Zero pivot encountered in Thomas algorithm");

    result.solution[n-1] = d[n-1] / b[n-1];
    result.flop_count += 1; // division

    for (int i = n-2; i >= 0; --i) {
        result.solution[i] = (d[i] - c[i] * result.solution[i+1]) / b[i];
        result.flop_count += 3; // multiplication, subtraction, division
    }

    result.iterations = 1;  // Direct method
    result.converged = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical

