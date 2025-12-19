#include "../../include/interpolation.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>

namespace numerical {

InterpolationResult lagrange_polynomial(
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
    result.interpolated_values.resize(x_query.size());
    
    // For each query point, compute Lagrange interpolation
    for (size_t q = 0; q < x_query.size(); ++q) {
        double x = x_query[q];
        double sum = 0.0;
        
        // Lagrange formula: L(x) = sum(y_i * L_i(x))
        // where L_i(x) = product((x - x_j) / (x_i - x_j)) for j != i
        for (size_t i = 0; i < n; ++i) {
            double L_i = 1.0;
            
            // Compute L_i(x) = product((x - x_j) / (x_i - x_j))
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    L_i *= (x - x_data[j]) / (x_data[i] - x_data[j]);
                    result.flop_count += 3; // subtract, subtract, divide, multiply
                }
            }
            
            sum += y_data[i] * L_i;
            result.flop_count += 1; // multiply, add
        }
        
        result.interpolated_values[q] = sum;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    result.success = true;
    
    return result;
}

} // namespace numerical

