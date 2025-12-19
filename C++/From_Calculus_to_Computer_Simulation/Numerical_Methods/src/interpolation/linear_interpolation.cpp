#include "../../include/interpolation.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace numerical {

InterpolationResult linear_interpolation(
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
    
    // Check if x_data is sorted
    for (size_t i = 1; i < x_data.size(); ++i) {
        if (x_data[i] <= x_data[i-1]) {
            throw std::invalid_argument("x_data must be strictly increasing");
        }
    }
    
    result.interpolated_values.resize(x_query.size());
    
    // Perform linear interpolation for each query point
    for (size_t q = 0; q < x_query.size(); ++q) {
        double x = x_query[q];
        
        // Handle extrapolation
        if (x < x_data[0]) {
            // Extrapolate to the left using first two points
            double slope = (y_data[1] - y_data[0]) / (x_data[1] - x_data[0]);
            result.interpolated_values[q] = y_data[0] + slope * (x - x_data[0]);
            result.flop_count += 4; // subtract, divide, subtract, add
        } else if (x > x_data.back()) {
            // Extrapolate to the right using last two points
            size_t n = x_data.size();
            double slope = (y_data[n-1] - y_data[n-2]) / (x_data[n-1] - x_data[n-2]);
            result.interpolated_values[q] = y_data[n-1] + slope * (x - x_data[n-1]);
            result.flop_count += 4;
        } else {
            // Find the interval containing x
            size_t i = 0;
            for (size_t j = 0; j < x_data.size() - 1; ++j) {
                if (x >= x_data[j] && x <= x_data[j+1]) {
                    i = j;
                    break;
                }
            }
            
            // Linear interpolation: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            double x0 = x_data[i];
            double x1 = x_data[i+1];
            double y0 = y_data[i];
            double y1 = y_data[i+1];
            
            double t = (x - x0) / (x1 - x0);
            result.interpolated_values[q] = y0 + t * (y1 - y0);
            result.flop_count += 5; // subtract, subtract, divide, subtract, multiply, add
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    result.success = true;
    
    return result;
}

} // namespace numerical

