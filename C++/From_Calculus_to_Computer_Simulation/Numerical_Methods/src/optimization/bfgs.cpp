#include "../../include/optimization.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

OptimizationResult bfgs(
    ObjectiveFunc f,
    GradientFunc grad_f,
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
    std::vector<double> x = x0;
    
    // Initialize inverse Hessian approximation as identity matrix
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        H[i][i] = 1.0;
    }
    
    std::vector<double> grad_old = grad_f(x);
    result.flop_count += n;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute search direction: d = -H * grad
        std::vector<double> d(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                d[i] -= H[i][j] * grad_old[j];
                result.flop_count += 2;
            }
        }
        
        // Line search with backtracking
        double alpha = 1.0;
        double f_old = f(x);
        std::vector<double> x_new(n);
        double c = 0.5;  // Armijo constant
        double tau = 0.5;  // Backtracking factor
        
        // Compute dot product for Armijo condition
        double dot_product = 0.0;
        for (size_t i = 0; i < n; ++i) {
            dot_product += grad_old[i] * d[i];
        }
        result.flop_count += n;
        
        // Backtracking line search
        for (int ls_iter = 0; ls_iter < 20; ++ls_iter) {
            for (size_t i = 0; i < n; ++i) {
                x_new[i] = x[i] + alpha * d[i];
            }
            double f_new = f(x_new);
            result.flop_count += 1;
            
            // Armijo condition: f(x + alpha*d) <= f(x) + c*alpha*grad^T*d
            if (f_new <= f_old + c * alpha * dot_product) {
                break;
            }
            alpha *= tau;
            result.flop_count += 1;
        }
        
        result.flop_count += n * 2; // n additions and n multiplications
        
        std::vector<double> grad_new = grad_f(x_new);
        result.flop_count += n;
        
        // Check convergence
        double grad_norm = 0.0;
        for (size_t i = 0; i < n; ++i) {
            grad_norm += grad_new[i] * grad_new[i];
            result.flop_count += 1;
        }
        grad_norm = std::sqrt(grad_norm);
        result.flop_count += 1;
        
        // Also check if function value changed significantly
        double f_new = f(x_new);
        double f_change = std::fabs(f_new - f_old);
        result.flop_count += 1;
        
        if (grad_norm < tolerance && f_change < tolerance) {
            result.converged = true;
            x = x_new;
            result.optimal_point = x;
            result.optimal_value = f_new;
            result.iterations = iter + 1;
            break;
        }
        
        // BFGS update
        std::vector<double> s(n), y(n);
        for (size_t i = 0; i < n; ++i) {
            s[i] = x_new[i] - x[i];
            y[i] = grad_new[i] - grad_old[i];
            result.flop_count += 2;
        }
        
        double rho = 0.0;
        for (size_t i = 0; i < n; ++i) {
            rho += s[i] * y[i];
            result.flop_count += 1;
        }
        
        if (std::fabs(rho) > 1e-10) {
            rho = 1.0 / rho;
            
            // Update H using BFGS formula: H_new = H + rho*(s*s^T - H*y*s^T - s*y^T*H)
            // Simplified version
            std::vector<std::vector<double>> H_new = H;
            std::vector<double> Hy(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    Hy[i] += H[i][j] * y[j];
                    result.flop_count += 1;
                }
            }
            
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    H_new[i][j] += rho * (s[i] * s[j] - Hy[i] * s[j] - s[i] * Hy[j]);
                    result.flop_count += 4;
                }
            }
            H = H_new;
        }
        
        x = x_new;
        grad_old = grad_new;
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

