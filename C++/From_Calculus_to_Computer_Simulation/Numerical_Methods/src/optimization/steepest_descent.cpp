#include "../../include/optimization.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

OptimizationResult steepest_descent(
    ObjectiveFunc f,
    GradientFunc grad_f,
    const std::vector<double>& x0,
    double step_size,
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
    
    std::vector<double> x = x0;
    size_t n = x.size();
    
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> grad = grad_f(x);
        result.flop_count += n;
        
        // Compute gradient norm
        double grad_norm = 0.0;
        for (size_t i = 0; i < n; ++i) {
            grad_norm += grad[i] * grad[i];
            result.flop_count += 1;
        }
        grad_norm = std::sqrt(grad_norm);
        result.flop_count += 1;
        
        // Check convergence
        if (grad_norm < tolerance) {
            result.converged = true;
            break;
        }
        
        // Update: x = x - step_size * grad
        for (size_t i = 0; i < n; ++i) {
            x[i] = x[i] - step_size * grad[i];
            result.flop_count += 2;
        }
        
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

