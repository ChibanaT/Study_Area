#include "../../include/optimization.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

OptimizationResult sqp(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<ConstraintFunc>& constraints,
    const std::vector<GradientFunc>& grad_constraints,
    const std::vector<double>& x0,
    double tolerance,
    int max_iter
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    OptimizationResult result;
    result.converged = false;
    result.flop_count = 0;
    result.iterations = 0;
    result.optimal_point = x0;
    result.optimal_value = f(x0);
    result.flop_count += 1;
    
    size_t n = x0.size();
    size_t m = constraints.size();
    std::vector<double> x = x0;
    
    // Sequential Quadratic Programming
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> grad = grad_f(x);
        result.flop_count += n;
        
        // Check constraint violations
        std::vector<double> constraint_values(m);
        for (size_t i = 0; i < m; ++i) {
            constraint_values[i] = constraints[i](x);
            result.flop_count += 1;
        }
        
        // Compute constraint gradients
        std::vector<std::vector<double>> constraint_grads(m);
        for (size_t i = 0; i < m; ++i) {
            constraint_grads[i] = grad_constraints[i](x);
            result.flop_count += n;
        }
        
        // Solve QP subproblem (simplified)
        // Minimize: 0.5 * d^T * H * d + grad^T * d
        // Subject to: constraint_grads[i]^T * d + constraint_values[i] <= 0
        
        // Use gradient projection method
        std::vector<double> d(n);
        for (size_t i = 0; i < n; ++i) {
            d[i] = -grad[i];
            result.flop_count += 1;
        }
        
        // Project onto feasible set
        for (size_t i = 0; i < m; ++i) {
            if (constraint_values[i] > -tolerance) {
                // Project d onto constraint boundary
                double norm = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    norm += constraint_grads[i][j] * constraint_grads[i][j];
                    result.flop_count += 1;
                }
                if (norm > 1e-10) {
                    double proj = 0.0;
                    for (size_t j = 0; j < n; ++j) {
                        proj += d[j] * constraint_grads[i][j];
                        result.flop_count += 1;
                    }
                    for (size_t j = 0; j < n; ++j) {
                        d[j] -= (proj / norm) * constraint_grads[i][j];
                        result.flop_count += 2;
                    }
                }
            }
        }
        
        // Line search
        double alpha = 1.0;
        std::vector<double> x_new(n);
        for (size_t i = 0; i < n; ++i) {
            x_new[i] = x[i] + alpha * d[i];
            result.flop_count += 2;
        }
        
        // Check convergence
        double grad_norm = 0.0;
        for (size_t i = 0; i < n; ++i) {
            grad_norm += grad[i] * grad[i];
            result.flop_count += 1;
        }
        grad_norm = std::sqrt(grad_norm);
        result.flop_count += 1;
        
        double max_violation = 0.0;
        for (size_t i = 0; i < m; ++i) {
            max_violation = std::max(max_violation, constraint_values[i]);
        }
        
        if (grad_norm < tolerance && max_violation < tolerance) {
            result.converged = true;
            x = x_new;
            result.optimal_point = x;
            result.optimal_value = f(x);
            result.iterations = iter + 1;
            break;
        }
        
        x = x_new;
        result.optimal_point = x;
        result.optimal_value = f(x);
        result.iterations = iter + 1;
        result.history.push_back(x);
        result.flop_count += 1;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

