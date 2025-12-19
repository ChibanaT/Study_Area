#pragma once
#include <functional>
#include <vector>

namespace numerical {

//----------------------------
// Struct to hold ODE solver results
//----------------------------
struct ODEResult {
    std::vector<double> time_points;          // Time points where solution was computed
    std::vector<std::vector<double>> solution; // Solution values at each time point (for systems)
    std::vector<double> solution_single;      // Solution values for single ODE
    int steps;                                 // Number of steps taken
    double cpu_time_sec;                       // CPU time in seconds
    long long flop_count;                      // Floating-point operation count
    bool converged;                            // Convergence flag
    double max_error;                          // Maximum error estimate (for adaptive methods)
};

//----------------------------
// Function types for ODEs
//----------------------------
using ODEFunc = std::function<double(double, double)>;  // dy/dt = f(t, y)
using ODESystemFunc = std::function<std::vector<double>(double, const std::vector<double>&)>;  // For systems

//----------------------------
// ODE solver method declarations
//----------------------------

// Euler method: First-order explicit method
ODEResult euler(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h
);

// Runge-Kutta 4th order: Classic RK4 method
ODEResult rk4(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h
);

// Runge-Kutta-Fehlberg 4(5): Adaptive step-size method
ODEResult rk45(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h_initial,
    double tolerance = 1e-6
);

// Runge-Kutta 7(8): High-order adaptive method
ODEResult rk78(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h_initial,
    double tolerance = 1e-6
);

// Backward Differentiation Formula (BDF): Implicit method for stiff ODEs
ODEResult bdf(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h,
    int order = 2
);

// Radau method: Implicit Runge-Kutta method for stiff systems
ODEResult radau(
    ODEFunc f,
    double t0,
    double y0,
    double t_end,
    double h,
    double tolerance = 1e-6
);

// System versions (for ODE systems)
ODEResult euler_system(
    ODESystemFunc f,
    double t0,
    const std::vector<double>& y0,
    double t_end,
    double h
);

ODEResult rk4_system(
    ODESystemFunc f,
    double t0,
    const std::vector<double>& y0,
    double t_end,
    double h
);

} // namespace numerical

