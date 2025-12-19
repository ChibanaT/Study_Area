#include "../../include/optimization.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

namespace numerical {

OptimizationResult genetic_algorithm(
    ObjectiveFunc f,
    int dimension,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int population_size,
    int max_generations,
    double mutation_rate,
    double crossover_rate
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
    
    // Initialize population
    std::vector<std::vector<double>> population(population_size);
    std::vector<double> fitness(population_size);
    
    for (int i = 0; i < population_size; ++i) {
        population[i].resize(dimension);
        for (int j = 0; j < dimension; ++j) {
            std::uniform_real_distribution<double> dist(lower_bounds[j], upper_bounds[j]);
            population[i][j] = dist(gen);
        }
        fitness[i] = f(population[i]);
        result.flop_count += 1;
    }
    
    // Find best individual
    int best_idx = 0;
    for (int i = 1; i < population_size; ++i) {
        if (fitness[i] < fitness[best_idx]) {
            best_idx = i;
        }
    }
    
    result.optimal_point = population[best_idx];
    result.optimal_value = fitness[best_idx];
    
    // Evolution loop
    for (int generation = 0; generation < max_generations; ++generation) {
        // Selection (tournament selection)
        std::vector<std::vector<double>> new_population(population_size);
        std::vector<double> new_fitness(population_size);
        
        for (int i = 0; i < population_size; ++i) {
            // Tournament selection
            std::uniform_int_distribution<int> dist_idx(0, population_size - 1);
            int idx1 = dist_idx(gen);
            int idx2 = dist_idx(gen);
            int parent1 = (fitness[idx1] < fitness[idx2]) ? idx1 : idx2;
            
            idx1 = dist_idx(gen);
            idx2 = dist_idx(gen);
            int parent2 = (fitness[idx1] < fitness[idx2]) ? idx1 : idx2;
            
            // Crossover
            std::uniform_real_distribution<double> dist_cross(0.0, 1.0);
            if (dist_cross(gen) < crossover_rate) {
                new_population[i] = population[parent1];
                int crossover_point = dist_idx(gen) % dimension;
                for (int j = crossover_point; j < dimension; ++j) {
                    new_population[i][j] = population[parent2][j];
                }
            } else {
                new_population[i] = population[parent1];
            }
            
            // Mutation
            std::uniform_real_distribution<double> dist_mut(0.0, 1.0);
            for (int j = 0; j < dimension; ++j) {
                if (dist_mut(gen) < mutation_rate) {
                    std::uniform_real_distribution<double> dist_val(lower_bounds[j], upper_bounds[j]);
                    new_population[i][j] = dist_val(gen);
                }
            }
            
            new_fitness[i] = f(new_population[i]);
            result.flop_count += 1;
        }
        
        population = new_population;
        fitness = new_fitness;
        
        // Update best
        for (int i = 0; i < population_size; ++i) {
            if (fitness[i] < result.optimal_value) {
                result.optimal_point = population[i];
                result.optimal_value = fitness[i];
            }
        }
        
        result.iterations = generation + 1;
        result.history.push_back(result.optimal_point);
    }
    
    result.converged = true;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    result.cpu_time_sec = duration.count() * 1e-9;
    
    return result;
}

} // namespace numerical

