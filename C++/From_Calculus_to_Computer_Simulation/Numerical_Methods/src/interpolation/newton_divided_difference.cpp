#include "../../include/interpolation.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>

namespace numerical {

InterpolationResult newton_divided_difference(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& x_query
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
    
    size_t n = x_data.size();
    
    // Compute divided differences table
    std::vector<std::vector<double>> div_diff(n, std::vector<double>(n));
    
    // Initialize first column with y values
    for (size_t i = 0; i < n; ++i) {
        div_diff[i][0] = y_data[i];
    }
    
    // Compute divided differences
    for (size_t j = 1; j < n; ++j) {
        for (size_t i = 0; i < n - j; ++i) {
            div_diff[i][j] = (div_diff[i+1][j-1] - div_diff[i][j-1]) / (x_data[i+j] - x_data[i]);
            result.flop_count += 3; // subtract, subtract, divide
        }
    }
    
    // Store coefficients (first row of divided differences table)
    result.coefficients.resize(n);
    for (size_t i = 0; i < n; ++i) {
        result.coefficients[i] = div_diff[0][i];
    }
    
    // Evaluate Newton polynomial for each query point
    result.interpolated_values.resize(x_query.size());
    
    for (size_t q = 0; q < x_query.size(); ++q) {
        double x = x_query[q];
        double value = div_diff[0][0];  // f[x0]
        double product = 1.0;
        
        // Newton form: N(x) = f[x0] + f[x0,x1](x-x0) + f[x0,x1,x2](x-x0)(x-x1) + ...
        for (size_t i = 1; i < n; ++i) {
            product *= (x - x_data[i-1]);
            value += div_diff[0][i] * product;
            result.flop_count += 3; // subtract, multiply, multiply, add
        }
        
        result.interpolated_values[q] = value;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    result.success = true;
    
    return result;
}

} // namespace numerical

