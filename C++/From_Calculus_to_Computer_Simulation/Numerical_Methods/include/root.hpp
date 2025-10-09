#pragma once
#include <functional>
#include <vector>

namespace numerical {

//----------------------------
// Struct to hold root-finding results
//----------------------------
struct RootResult {
    double root;                     // Last computed root
    int iterations;                  // Number of iterations
    bool converged;                  // Convergence flag
    std::vector<double> history;     // History of approximations
    double cpu_time_sec;             // CPU time in seconds
    long long flop_count;            // Count of floating-point operations
};

//----------------------------
// Function type
//----------------------------
using Func = std::function<double(double)>;

//----------------------------
// Root-finding method declarations
//----------------------------
RootResult bisect(const Func &f, double a, double b, double tol, int max_iter);
RootResult false_position(const Func &f, double a, double b, double tol, int max_iter);
RootResult modified_false_position(const Func &f, double a, double b, double tol, int max_iter);
RootResult secant(const Func &f, double x0, double x1, double tol, int max_iter);
RootResult newton_raphson(const Func &f, const Func &df, double x0, double tol, int max_iter);
RootResult modified_newton(const Func &f, const Func &df, double x0, double tol, int max_iter);
RootResult halley(const Func &f, const Func &df, const Func &d2f, double x0, double tol, int max_iter);
RootResult brent(const Func &f, double a, double b, double tol, int max_iter);
RootResult bairstow(const std::vector<double>& coeffs, double r, double s, double tol, int max_iter);
RootResult muller(const Func &f, double x0, double x1, double x2, double tol, int max_iter);

//----------------------------
// Optional helper functions
//----------------------------
double horner(const std::vector<double>& coeffs, double x);          // Evaluate polynomial
double derivative(const std::vector<double>& coeffs, double x);     // First derivative
double second_derivative(const std::vector<double>& coeffs, double x); // Second derivative

} // namespace numerical
