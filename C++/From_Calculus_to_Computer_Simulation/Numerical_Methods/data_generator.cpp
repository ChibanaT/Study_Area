#include <iostream>
#include <fstream>
#include <cmath>   
#include <cstdlib>
#include <vector>
#include <ctime> 

// Function to generate a normally distributed random number (Box-Muller transform approximation)
double generateNormalNoise(double mean, double stddev) {
    double u1 = (double)std::rand() / RAND_MAX;
    double u2 = (double)std::rand() / RAND_MAX;
    double rand_std_normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    return mean + stddev * rand_std_normal;
}

// Simple Bubble Sort implementation (fewer libraries)
void sortVector(std::vector<double>& vec) {
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (vec[j] > vec[j + 1]) {
                double temp = vec[j];
                vec[j] = vec[j + 1];
                vec[j + 1] = temp;
            }
        }
    }
}

void generateSparseNoisyData(int num_points, const std::string& filename) {
    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::vector<double> x_values(num_points);
    std::vector<double> y_noisy_values(num_points);

    // Generate X values (using uniform random and sorting to simulate desired distribution)
    for (int i = 0; i < num_points; ++i) {
        // Generate uniform random and add offset
        x_values[i] = (double)std::rand() / RAND_MAX * 25.0 + 0.1; 
    }

    // Sort X values manually to simulate increasing sparsity (densidade variável)
    sortVector(x_values);

    // Generate Y values with manual normal noise
    for (int i = 0; i < num_points; ++i) {
        double true_y = std::log(x_values[i]);
        double noise = generateNormalNoise(0.0, 0.2); // mean=0.0, stddev=0.2
        y_noisy_values[i] = true_y + noise;
    }

    // Write data to the file
    std::ofstream output_file(filename);

    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    

    for (int i = 0; i < num_points; ++i) {
        
        output_file << x_values[i] << ";" << y_noisy_values[i] << "\n";
    }

    output_file.close();
    std::cout << "Successfully generated " << num_points << " data points in " << filename << std::endl;
}

int main() {
    generateSparseNoisyData(500, "data.dat");
    return 0;
}
