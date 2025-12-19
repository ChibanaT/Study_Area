#include "../../include/ode_solver.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>

namespace numerical {

ODEResult rk4(
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
    
    double t = t0;
    double y = y0;
    result.time_points.push_back(t);
    result.solution_single.push_back(y);
    
    // Runge-Kutta 4th order method
    while (t < t_end) {
        double h_actual = h;
        if (t + h > t_end) {
            h_actual = t_end - t;
        }
        
        double k1 = f(t, y);
        double k2 = f(t + h_actual/2.0, y + h_actual*k1/2.0);
        double k3 = f(t + h_actual/2.0, y + h_actual*k2/2.0);
        double k4 = f(t + h_actual, y + h_actual*k3);
        result.flop_count += 4; // 4 function evaluations
        
        y = y + (h_actual/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t = t + h_actual;
        result.flop_count += 10; // Various arithmetic operations
        
        result.time_points.push_back(t);
        result.solution_single.push_back(y);
        result.steps++;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

ODEResult rk4_system(
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
        std::vector<double> y_temp(n);
        
        for (size_t i = 0; i < n; ++i) y_temp[i] = y[i] + h*k1[i]/2.0;
        std::vector<double> k2 = f(t + h/2.0, y_temp);
        
        for (size_t i = 0; i < n; ++i) y_temp[i] = y[i] + h*k2[i]/2.0;
        std::vector<double> k3 = f(t + h/2.0, y_temp);
        
        for (size_t i = 0; i < n; ++i) y_temp[i] = y[i] + h*k3[i];
        std::vector<double> k4 = f(t + h, y_temp);
        
        result.flop_count += 4 * n; // 4 function evaluations
        
        for (size_t i = 0; i < n; ++i) {
            y[i] = y[i] + (h/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
            result.flop_count += 6;
        }
        
        t = t + h;
        result.flop_count += 1;
        
        result.time_points.push_back(t);
        result.solution.push_back(y);
        result.steps++;
        
        if (t > t_end) {
            double h_last = t_end - (t - h);
            // Recompute last step with correct step size
            k1 = f(t - h, result.solution[result.solution.size() - 2]);
            for (size_t i = 0; i < n; ++i) y_temp[i] = result.solution[result.solution.size() - 2][i] + h_last*k1[i]/2.0;
            k2 = f(t - h + h_last/2.0, y_temp);
            for (size_t i = 0; i < n; ++i) y_temp[i] = result.solution[result.solution.size() - 2][i] + h_last*k2[i]/2.0;
            k3 = f(t - h + h_last/2.0, y_temp);
            for (size_t i = 0; i < n; ++i) y_temp[i] = result.solution[result.solution.size() - 2][i] + h_last*k3[i];
            k4 = f(t_end, y_temp);
            for (size_t i = 0; i < n; ++i) {
                y[i] = result.solution[result.solution.size() - 2][i] + (h_last/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
            }
            t = t_end;
            result.time_points.back() = t;
            result.solution.back() = y;
            result.flop_count += 14 * n;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

