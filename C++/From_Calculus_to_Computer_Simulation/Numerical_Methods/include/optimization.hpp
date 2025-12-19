#pragma once
#include <functional>
#include <vector>
#include <string>

namespace numerical {

//----------------------------
// Struct to hold optimization results
//----------------------------
struct OptimizationResult {
    std::vector<double> optimal_point;        // Optimal point found
    double optimal_value;                     // Optimal function value
    int iterations;                           // Number of iterations
    bool converged;                            // Convergence flag
    double cpu_time_sec;                      // CPU time in seconds
    long long flop_count;                     // Floating-point operation count
    std::vector<std::vector<double>> history; // History of iterations
};

//----------------------------
// Function types for optimization
//----------------------------
using ObjectiveFunc = std::function<double(const std::vector<double>&)>;
using GradientFunc = std::function<std::vector<double>(const std::vector<double>&)>;
using ConstraintFunc = std::function<double(const std::vector<double>&)>;  // g(x) <= 0

//----------------------------
// Optimization method declarations
//----------------------------

// Steepest Descent: First-order gradient method
OptimizationResult steepest_descent(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<double>& x0,
    double step_size,
    double tolerance,
    int max_iter
);

// BFGS: Quasi-Newton method with BFGS update
OptimizationResult bfgs(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<double>& x0,
    double tolerance,
    int max_iter
);

// Simplex method: Linear programming solver
OptimizationResult simplex(
    const std::vector<double>& c,                    // Objective function coefficients
    const std::vector<std::vector<double>>& A,      // Constraint matrix
    const std::vector<double>& b,                   // Constraint bounds
    const std::vector<std::string>& constraint_type, // "leq", "eq", "geq"
    const std::vector<double>& x0                    // Initial guess
);

// Sequential Quadratic Programming (SQP): For constrained optimization
OptimizationResult sqp(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<ConstraintFunc>& constraints,  // g_i(x) <= 0
    const std::vector<GradientFunc>& grad_constraints,
    const std::vector<double>& x0,
    double tolerance,
    int max_iter
);

// Lagrange Multipliers: For equality constrained optimization
OptimizationResult lagrange_multipliers(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<ConstraintFunc>& equality_constraints,  // h_i(x) = 0
    const std::vector<GradientFunc>& grad_equality_constraints,
    const std::vector<double>& x0,
    double tolerance,
    int max_iter
);

// Genetic Algorithm: Population-based metaheuristic
OptimizationResult genetic_algorithm(
    ObjectiveFunc f,
    int dimension,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int population_size = 50,
    int max_generations = 100,
    double mutation_rate = 0.1,
    double crossover_rate = 0.8
);

// Particle Swarm Optimization (PSO): Swarm intelligence method
OptimizationResult pso(
    ObjectiveFunc f,
    int dimension,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int swarm_size = 30,
    int max_iterations = 100,
    double inertia_weight = 0.7,
    double cognitive_coeff = 1.5,
    double social_coeff = 1.5
);

} // namespace numerical

