#include "../../include/ode_solver.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace numerical {

ODEResult rk45(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h_initial,
    double tolerance
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    ODEResult result;
    result.converged = true;
    result.flop_count = 0;
    result.steps = 0;
    result.max_error = 0.0;
    
    if (h_initial <= 0.0) {
        throw std::invalid_argument("Initial step size must be positive");
    }
    
    if (t_end <= t0) {
        throw std::invalid_argument("t_end must be greater than t0");
    }
    
    double t = t0;
    double y = y0;
    double h = h_initial;
    
    result.time_points.push_back(t);
    result.solution_single.push_back(y);
    
    // Runge-Kutta-Fehlberg 4(5) coefficients
    const double c2 = 1.0/4.0, c3 = 3.0/8.0, c4 = 12.0/13.0, c5 = 1.0, c6 = 1.0/2.0;
    const double a21 = 1.0/4.0;
    const double a31 = 3.0/32.0, a32 = 9.0/32.0;
    const double a41 = 1932.0/2197.0, a42 = -7200.0/2197.0, a43 = 7296.0/2197.0;
    const double a51 = 439.0/216.0, a52 = -8.0, a53 = 3680.0/513.0, a54 = -845.0/4104.0;
    const double a61 = -8.0/27.0, a62 = 2.0, a63 = -3544.0/2565.0, a64 = 1859.0/4104.0, a65 = -11.0/40.0;
    
    const double b41 = 25.0/216.0, b42 = 0.0, b43 = 1408.0/2565.0, b44 = 2197.0/4104.0, b45 = -1.0/5.0, b46 = 0.0;
    const double b51 = 16.0/135.0, b52 = 0.0, b53 = 6656.0/12825.0, b54 = 28561.0/56430.0, b55 = -9.0/50.0, b56 = 2.0/55.0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        double k1 = f(t, y);
        double k2 = f(t + c2*h, y + h*a21*k1);
        double k3 = f(t + c3*h, y + h*(a31*k1 + a32*k2));
        double k4 = f(t + c4*h, y + h*(a41*k1 + a42*k2 + a43*k3));
        double k5 = f(t + c5*h, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4));
        double k6 = f(t + c6*h, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5));
        
        result.flop_count += 6; // 6 function evaluations
        
        // 4th order estimate
        double y4 = y + h*(b41*k1 + b42*k2 + b43*k3 + b44*k4 + b45*k5 + b46*k6);
        // 5th order estimate
        double y5 = y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4 + b55*k5 + b56*k6);
        
        double error = std::fabs(y5 - y4);
        result.max_error = std::max(result.max_error, error);
        
        if (error < tolerance * h) {
            y = y5;  // Use higher order estimate
            t = t + h;
            result.time_points.push_back(t);
            result.solution_single.push_back(y);
            result.steps++;
            
            // Adapt step size
            if (error > 0.0) {
                h = 0.9 * h * std::pow(tolerance * h / error, 0.2);
            } else {
                h = 2.0 * h;
            }
            h = std::min(h, t_end - t);
            result.flop_count += 20;
        } else {
            // Reject step and reduce h
            h = 0.9 * h * std::pow(tolerance * h / error, 0.25);
            result.flop_count += 3;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

