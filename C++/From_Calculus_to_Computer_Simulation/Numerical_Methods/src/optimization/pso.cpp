#include "../../include/optimization.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <random>

namespace numerical {

OptimizationResult pso(
    ObjectiveFunc f,
    int dimension,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int swarm_size,
    int max_iterations,
    double inertia_weight,
    double cognitive_coeff,
    double social_coeff
) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    OptimizationResult result;
    result.converged = false;
    result.flop_count = 0;
    result.iterations = 0;
    
    if (lower_bounds.size() != dimension || upper_bounds.size() != dimension) {
        throw std::invalid_argument("Bounds must match dimension");
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Initialize particles
    std::vector<std::vector<double>> positions(swarm_size);
    std::vector<std::vector<double>> velocities(swarm_size);
    std::vector<std::vector<double>> best_positions(swarm_size);
    std::vector<double> best_values(swarm_size);
    
    std::vector<double> global_best_position(dimension);
    double global_best_value = 1e10;
    
    for (int i = 0; i < swarm_size; ++i) {
        positions[i].resize(dimension);
        velocities[i].resize(dimension);
        best_positions[i].resize(dimension);
        
        for (int j = 0; j < dimension; ++j) {
            std::uniform_real_distribution<double> dist_pos(lower_bounds[j], upper_bounds[j]);
            std::uniform_real_distribution<double> dist_vel(-1.0, 1.0);
            positions[i][j] = dist_pos(gen);
            velocities[i][j] = dist_vel(gen);
        }
        
        best_positions[i] = positions[i];
        best_values[i] = f(positions[i]);
        result.flop_count += 1;
        
        if (best_values[i] < global_best_value) {
            global_best_value = best_values[i];
            global_best_position = positions[i];
        }
    }
    
    result.optimal_point = global_best_position;
    result.optimal_value = global_best_value;
    
    // PSO main loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::uniform_real_distribution<double> dist_r(0.0, 1.0);
        
        for (int i = 0; i < swarm_size; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double r1 = dist_r(gen);
                double r2 = dist_r(gen);
                
                // Update velocity
                velocities[i][j] = inertia_weight * velocities[i][j] +
                    cognitive_coeff * r1 * (best_positions[i][j] - positions[i][j]) +
                    social_coeff * r2 * (global_best_position[j] - positions[i][j]);
                result.flop_count += 5;
                
                // Update position
                positions[i][j] += velocities[i][j];
                result.flop_count += 1;
                
                // Apply bounds
                positions[i][j] = std::max(lower_bounds[j], std::min(upper_bounds[j], positions[i][j]));
            }
            
            // Evaluate fitness
            double fitness = f(positions[i]);
            result.flop_count += 1;
            
            // Update personal best
            if (fitness < best_values[i]) {
                best_positions[i] = positions[i];
                best_values[i] = fitness;
                
                // Update global best
                if (fitness < global_best_value) {
                    global_best_value = fitness;
                    global_best_position = positions[i];
                }
            }
        }
        
        result.optimal_point = global_best_position;
        result.optimal_value = global_best_value;
        result.iterations = iter + 1;
        result.history.push_back(global_best_position);
    }
    
    result.converged = true;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

