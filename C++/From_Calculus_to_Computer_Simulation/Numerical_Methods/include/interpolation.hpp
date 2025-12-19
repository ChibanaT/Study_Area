#pragma once
#include <vector>

namespace numerical {

//----------------------------
// Struct to hold interpolation results
//----------------------------
struct InterpolationResult {
    std::vector<double> interpolated_values;  // Interpolated values
    std::vector<double> coefficients;         // Polynomial coefficients (if applicable)
    double cpu_time_sec;                      // CPU time in seconds
    long long flop_count;                     // Floating-point operation count
    bool success;                             // Success flag
};

//----------------------------
// Interpolation method declarations
//----------------------------

// Linear interpolation: Simple piecewise linear interpolation
InterpolationResult linear_interpolation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query
);

// Lagrange polynomial: Interpolates using Lagrange basis polynomials
InterpolationResult lagrange_polynomial(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query
);

// Newton divided difference: Interpolates using Newton's divided difference formula
InterpolationResult newton_divided_difference(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query
);

// Cubic spline: Smooth piecewise cubic interpolation with continuous derivatives
InterpolationResult cubic_spline(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query,
    double left_boundary_derivative = 0.0,    // Natural spline if not specified
    double right_boundary_derivative = 0.0    // Natural spline if not specified
);

} // namespace numerical

