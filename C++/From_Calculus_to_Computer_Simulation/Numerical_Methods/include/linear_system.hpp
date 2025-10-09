#pragma once
#include <vector>

namespace numerical {

//----------------------------
// Struct to hold linear system results
//----------------------------
struct LinearResult {
    std::vector<double> solution;   // Solution vector
    int iterations;                 // Number of iterations (for iterative methods)
    bool converged;                 // Convergence flag
    double cpu_time_sec;            // CPU time in seconds
    long long flop_count;           // Floating-point operation count
};

//----------------------------
// Linear system solver declarations
//----------------------------

// Direct methods
LinearResult cramer(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
LinearResult gaussian_elimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
LinearResult pivot_gauss(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
LinearResult gauss_jordan(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
LinearResult lu_decomposition(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
// LinearResult cholesky_decomposition(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
// LinearResult thomas_algorithm(const std::vector<std::vector<double>>& A, const std::vector<double>& b); // Tridiagonal

// // Iterative methods
// LinearResult jacobi(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int max_iter, double tol);
// LinearResult gauss_seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int max_iter, double tol);
// LinearResult sor(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int max_iter, double tol, double omega);
// LinearResult conjugate_gradient(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int max_iter, double tol);

} // namespace numerical
