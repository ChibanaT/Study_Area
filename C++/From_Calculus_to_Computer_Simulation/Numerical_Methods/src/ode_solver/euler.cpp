#include "../../include/ode_solver.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>

namespace numerical {

ODEResult euler(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    ODEResult result;
    result.converged = true;
    result.flop_count = 0;
    result.steps = 0;
    result.max_error = 0.0;
    
    if (h <= 0.0) {
        throw std::invalid_argument("Step size h must be positive");
    }
    
    if (t_end <= t0) {
        throw std::invalid_argument("t_end must be greater than t0");
    }
    
    // Initialize
    double t = t0;
    double y = y0;
    result.time_points.push_back(t);
    result.solution_single.push_back(y);
    
    // Euler method: y_{n+1} = y_n + h * f(t_n, y_n)
    while (t < t_end) {
        double h_actual = h;
        if (t + h > t_end) {
            h_actual = t_end - t;
        }
        
        double k1 = f(t, y);
        result.flop_count += 1; // function evaluation
        
        y = y + h_actual * k1;
        t = t + h_actual;
        result.flop_count += 3; // multiply, add, add
        
        result.time_points.push_back(t);
        result.solution_single.push_back(y);
        result.steps++;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

ODEResult euler_system(
    ODESystemFunc f,
    double t0,
    const std::vector<double>& y0,
    double t_end,
    double h
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    ODEResult result;
    result.converged = true;
    result.flop_count = 0;
    result.steps = 0;
    result.max_error = 0.0;
    
    if (h <= 0.0) {
        throw std::invalid_argument("Step size h must be positive");
    }
    
    if (t_end <= t0) {
        throw std::invalid_argument("t_end must be greater than t0");
    }
    
    size_t n = y0.size();
    std::vector<double> y = y0;
    double t = t0;
    
    result.time_points.push_back(t);
    result.solution.push_back(y);
    
    while (t < t_end) {
        std::vector<double> k1 = f(t, y);
        result.flop_count += n; // function evaluation
        
        for (size_t i = 0; i < n; ++i) {
            y[i] = y[i] + h * k1[i];
            result.flop_count += 2; // multiply, add
        }
        
        t = t + h;
        result.flop_count += 1;
        
        result.time_points.push_back(t);
        result.solution.push_back(y);
        result.steps++;
        
        if (t > t_end) {
            double h_last = t_end - (t - h);
            for (size_t i = 0; i < n; ++i) {
                y[i] = y[i] - h * k1[i] + h_last * k1[i];
                result.flop_count += 3;
            }
            t = t_end;
            result.time_points.back() = t;
            result.solution.back() = y;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

