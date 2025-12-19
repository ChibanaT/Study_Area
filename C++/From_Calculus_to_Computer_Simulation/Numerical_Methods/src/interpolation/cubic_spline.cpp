#include "../../include/interpolation.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>

namespace numerical {

InterpolationResult cubic_spline(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query,
    double left_boundary_derivative,
    double right_boundary_derivative
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    InterpolationResult result;
    result.success = false;
    result.flop_count = 0;
    result.interpolated_values.clear();
    
    // Validate inputs
    if (x_data.size() != y_data.size() || x_data.size() < 2) {
        throw std::invalid_argument("x_data and y_data must have same size and at least 2 points");
    }
    
    if (x_query.empty()) {
        return result;
    }
    
    // Check if x_data is sorted
    for (size_t i = 1; i < x_data.size(); ++i) {
        if (x_data[i] <= x_data[i-1]) {
            throw std::invalid_argument("x_data must be strictly increasing");
        }
    }
    
    size_t n = x_data.size();
    
    // Check if natural spline (both derivatives are zero and not explicitly set)
    bool natural_spline = (left_boundary_derivative == 0.0 && right_boundary_derivative == 0.0);
    
    // Compute second derivatives (s''(x)) using tridiagonal system
    std::vector<double> h(n - 1);  // Differences between x points
    std::vector<double> alpha(n);   // Right-hand side of tridiagonal system
    
    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = x_data[i+1] - x_data[i];
        result.flop_count += 1; // subtract
    }
    
    // Build tridiagonal system for natural spline or clamped spline
    std::vector<double> a(n), b(n), c(n);
    
    if (natural_spline) {
        // Natural spline: s''(x0) = s''(xn) = 0
        b[0] = 1.0;
        c[0] = 0.0;
        alpha[0] = 0.0;
        
        for (size_t i = 1; i < n - 1; ++i) {
            a[i] = h[i-1];
            b[i] = 2.0 * (h[i-1] + h[i]);
            c[i] = h[i];
            alpha[i] = 3.0 * ((y_data[i+1] - y_data[i]) / h[i] - (y_data[i] - y_data[i-1]) / h[i-1]);
            result.flop_count += 8; // multiply, subtract, divide, subtract, divide, subtract, multiply
        }
        
        a[n-1] = 0.0;
        b[n-1] = 1.0;
        alpha[n-1] = 0.0;
    } else {
        // Clamped spline: s'(x0) and s'(xn) are specified
        b[0] = 2.0 * h[0];
        c[0] = h[0];
        alpha[0] = 3.0 * ((y_data[1] - y_data[0]) / h[0] - left_boundary_derivative);
        result.flop_count += 4;
        
        for (size_t i = 1; i < n - 1; ++i) {
            a[i] = h[i-1];
            b[i] = 2.0 * (h[i-1] + h[i]);
            c[i] = h[i];
            alpha[i] = 3.0 * ((y_data[i+1] - y_data[i]) / h[i] - (y_data[i] - y_data[i-1]) / h[i-1]);
            result.flop_count += 8;
        }
        
        a[n-1] = h[n-2];
        b[n-1] = 2.0 * h[n-2];
        alpha[n-1] = 3.0 * (right_boundary_derivative - (y_data[n-1] - y_data[n-2]) / h[n-2]);
        result.flop_count += 4;
    }
    
    // Solve tridiagonal system using Thomas algorithm
    std::vector<double> c_prime(n);
    std::vector<double> d_prime(n);
    std::vector<double> second_deriv(n);
    
    c_prime[0] = c[0] / b[0];
    d_prime[0] = alpha[0] / b[0];
    result.flop_count += 2;
    
    for (size_t i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * c_prime[i-1];
        c_prime[i] = c[i] / denom;
        d_prime[i] = (alpha[i] - a[i] * d_prime[i-1]) / denom;
        result.flop_count += 6;
    }
    
    second_deriv[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; --i) {
        second_deriv[i] = d_prime[i] - c_prime[i] * second_deriv[i+1];
        result.flop_count += 2;
    }
    
    // Evaluate spline for each query point
    result.interpolated_values.resize(x_query.size());
    
    for (size_t q = 0; q < x_query.size(); ++q) {
        double x = x_query[q];
        
        // Find the interval
        size_t i = 0;
        if (x < x_data[0]) {
            i = 0;
        } else if (x >= x_data[n-1]) {
            i = n - 2;
        } else {
            for (size_t j = 0; j < n - 1; ++j) {
                if (x >= x_data[j] && x <= x_data[j+1]) {
                    i = j;
                    break;
                }
            }
        }
        
        // Evaluate cubic spline: s(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
        // Standard cubic spline formula
        double dx = x - x_data[i];
        double dx2 = dx * dx;
        double dx3 = dx2 * dx;
        
        // Coefficients for the cubic polynomial in the interval [x_i, x_{i+1}]
        // a_i = y_i (already have)
        // b_i = (y_{i+1} - y_i)/h_i - h_i*(2*m_i + m_{i+1})/3
        // c_i = m_i/2
        // d_i = (m_{i+1} - m_i)/(6*h_i)
        // where m_i = second_deriv[i]
        double a_coeff = y_data[i];
        double b_coeff = (y_data[i+1] - y_data[i]) / h[i] - h[i] * (2.0 * second_deriv[i] + second_deriv[i+1]) / 3.0;
        double c_coeff = second_deriv[i] / 2.0;
        double d_coeff = (second_deriv[i+1] - second_deriv[i]) / (6.0 * h[i]);
        
        result.interpolated_values[q] = a_coeff + b_coeff * dx + c_coeff * dx2 + d_coeff * dx3;
        result.flop_count += 15; // Various operations
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    result.success = true;
    
    return result;
}

} // namespace numerical

