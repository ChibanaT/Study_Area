#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <cstdio>   // popen, fprintf, fflush
#include <thread>   // sleep_for
#include <chrono>   // milliseconds
#include <algorithm>

struct Particle {
    double x, y, z;
    double vx, vy, vz;
    int source; // 0 = tube1, 1 = tube2
};

const double pi = 3.141592653589793;
const double dt = 0.01;          // time step (s)
const int total_steps = 1000;    // simulation steps
const double mu_air = 0.000018;  // air viscosity (kg/m*s)

// Tube parameters (converted to meters for physics, displayed in mm)
const double tube1_diam = 0.02;
const double tube1_angle = 30.0;
const double tube1_x = 0.0, tube1_y = 0.0;

const double tube2_diam = 0.02;
const double tube2_angle = 60.0;
const double tube2_x = 0.2, tube2_y = 0.0;  // 200mm

// Tank parameters (meters)
const double tank_diam = 0.4;       // 400mm
const double tank_height = 0.5;     // 500mm

// Brownian diffusion
const double sigma_diffusion = 0.0001; // m/s

// Random number generator
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> normal_dist(0.0, 1.0);
std::uniform_real_distribution<> uniform_dist(0.0, 1.0);

double randn() {
    return normal_dist(gen);
}

// 3D velocity field interpolation
struct Vec3 { double vx, vy, vz; };

Vec3 velocityField(double x, double y, double z) {
    double base_velocity = 0.5; // m/s
    
    double vx1 = cos(tube1_angle * pi / 180.0) * base_velocity;
    double vy1 = sin(tube1_angle * pi / 180.0) * base_velocity;
    double vx2 = cos(tube2_angle * pi / 180.0) * base_velocity;
    double vy2 = sin(tube2_angle * pi / 180.0) * base_velocity;

    double d1 = std::sqrt((x - tube1_x) * (x - tube1_x) + (y - tube1_y) * (y - tube1_y));
    double d2 = std::sqrt((x - tube2_x) * (x - tube2_x) + (y - tube2_y) * (y - tube2_y));
    
    double eps = 0.001;
    double w1 = 1.0 / (d1 + eps);
    double w2 = 1.0 / (d2 + eps);
    double wsum = w1 + w2;
    w1 /= wsum;
    w2 /= wsum;

    Vec3 v;
    v.vx = w1 * vx1 + w2 * vx2;
    v.vy = w1 * vy1 + w2 * vy2;
    v.vz = 0.0;
    
    // Add mixing circulation
    double cx = tank_diam / 2.0;
    double cy = tank_diam / 2.0;
    double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
    if (r > 0.01) {
        double circulation_strength = 0.1;
        v.vx += circulation_strength * (-(y - cy)) / r;
        v.vy += circulation_strength * (x - cx) / r;
    }
    
    return v;
}

// Collision / repulsion in 3D
void handleCollisions(std::vector<Particle> &particles, double particle_radius) {
    size_t N = particles.size();
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;
            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            
            if (dist < 2 * particle_radius && dist > 0) {
                double overlap = 0.5 * (2 * particle_radius - dist) / dist;
                particles[i].x -= dx * overlap;
                particles[i].y -= dy * overlap;
                particles[i].z -= dz * overlap;
                particles[j].x += dx * overlap;
                particles[j].y += dy * overlap;
                particles[j].z += dz * overlap;
            }
        }
    }
}

// Create tank and pipe geometry files
void createGeometry() {
    // Create tank walls (transparent box)
    std::ofstream tank_file("tank.dat");
    if (tank_file.is_open()) {
        double w = 400.0, h = 500.0; // mm
        
        // Bottom face
        tank_file << "0 0 0\n" << w << " 0 0\n" << w << " " << w << " 0\n" << "0 " << w << " 0\n" << "0 0 0\n\n";
        // Top face  
        tank_file << "0 0 " << h << "\n" << w << " 0 " << h << "\n" << w << " " << w << " " << h << "\n" << "0 " << w << " " << h << "\n" << "0 0 " << h << "\n\n";
        // Vertical edges
        tank_file << "0 0 0\n" << "0 0 " << h << "\n\n";
        tank_file << w << " 0 0\n" << w << " 0 " << h << "\n\n";
        tank_file << w << " " << w << " 0\n" << w << " " << w << " " << h << "\n\n";
        tank_file << "0 " << w << " 0\n" << "0 " << w << " " << h << "\n\n";
        tank_file.close();
    }
    
    // Create pipe 1 (tube1)
    std::ofstream pipe1_file("pipe1.dat");
    if (pipe1_file.is_open()) {
        double x_start = tube1_x * 1000; // Convert to mm
        double y_start = tube1_y * 1000;
        double radius = tube1_diam * 1000 / 2.0; // 10mm radius
        double length = 80.0; // mm pipe length
        
        // Pipe direction
        double dx = cos(tube1_angle * pi / 180.0);
        double dy = sin(tube1_angle * pi / 180.0);
        
        // Create cylindrical pipe segments
        for (int i = 0; i <= 20; i++) {
            double angle = i * 2 * pi / 20;
            double x_offset = radius * cos(angle);
            double y_offset = radius * sin(angle);
            
            // Rotate offset by pipe angle
            double x_rot = x_offset * cos(tube1_angle * pi / 180.0) - y_offset * sin(tube1_angle * pi / 180.0);
            double y_rot = x_offset * sin(tube1_angle * pi / 180.0) + y_offset * cos(tube1_angle * pi / 180.0);
            
            // Start of pipe
            pipe1_file << (x_start - length * dx + x_rot) << " " << (y_start - length * dy + y_rot) << " 0\n";
            // End of pipe (at tank)
            pipe1_file << (x_start + x_rot) << " " << (y_start + y_rot) << " 0\n\n";
        }
        pipe1_file.close();
    }
    
    // Create pipe 2 (tube2)
    std::ofstream pipe2_file("pipe2.dat");
    if (pipe2_file.is_open()) {
        double x_start = tube2_x * 1000; // Convert to mm
        double y_start = tube2_y * 1000;
        double radius = tube2_diam * 1000 / 2.0; // 10mm radius
        double length = 80.0; // mm pipe length
        
        // Pipe direction
        double dx = cos(tube2_angle * pi / 180.0);
        double dy = sin(tube2_angle * pi / 180.0);
        
        // Create cylindrical pipe segments
        for (int i = 0; i <= 20; i++) {
            double angle = i * 2 * pi / 20;
            double x_offset = radius * cos(angle);
            double y_offset = radius * sin(angle);
            
            // Rotate offset by pipe angle
            double x_rot = x_offset * cos(tube2_angle * pi / 180.0) - y_offset * sin(tube2_angle * pi / 180.0);
            double y_rot = x_offset * sin(tube2_angle * pi / 180.0) + y_offset * cos(tube2_angle * pi / 180.0);
            
            // Start of pipe
            pipe2_file << (x_start - length * dx + x_rot) << " " << (y_start - length * dy + y_rot) << " 0\n";
            // End of pipe (at tank)
            pipe2_file << (x_start + x_rot) << " " << (y_start + y_rot) << " 0\n\n";
        }
        pipe2_file.close();
    }
    
    std::cout << "Created geometry files: tank.dat, pipe1.dat, pipe2.dat\n";
}

// Export particles to single file (coordinates in mm for visualization)
void exportParticles(const std::vector<Particle> &particles, const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open " << filename << " for writing\n";
        return;
    }
    
    for (const auto &p : particles) {
        // Convert to mm for visualization
        file << p.x * 1000 << " " << p.y * 1000 << " " << p.z * 1000 << " " << p.source << "\n";
    }
    file.close();
}

int main() {
    // Input: tube 1
    double flow_rate1, density1;
    std::cout << "Tube 1 particle flow rate (g/s): ";
    std::cin >> flow_rate1;
    std::cout << "Tube 1 particle density (g/cm^3): ";
    std::cin >> density1;

    // Input: tube 2
    double flow_rate2, density2;
    std::cout << "Tube 2 particle flow rate (g/s): ";
    std::cin >> flow_rate2;
    std::cout << "Tube 2 particle density (g/cm^3): ";
    std::cin >> density2;

    // Particle parameters
    double particle_radius = 0.001; // 1mm -> 0.001m
    
    // Convert densities to kg/m^3
    double density1_si = density1 * 1000;
    double density2_si = density2 * 1000;
    
    // Calculate particle masses in kg
    double particle_volume = (4.0/3.0) * pi * pow(particle_radius, 3);
    double particle_mass1 = particle_volume * density1_si;
    double particle_mass2 = particle_volume * density2_si;

    // Initialize particles with good separation
    std::vector<Particle> particles;
    int N_particles = 50; // Reduced for better visibility
    
    for (int i = 0; i < N_particles; i++) {
        double spread = 0.02; // 20mm spread
        
        // Tube 1 particles (red) - start near tube1
        particles.push_back({
            tube1_x + spread * uniform_dist(gen),
            tube1_y + spread * uniform_dist(gen),
            0.05 + 0.1 * uniform_dist(gen), // 50-150mm height
            0.1 + 0.05 * randn(), // Initial velocity toward center
            0.1 + 0.05 * randn(),
            0.02 * randn(),
            0
        });
        
        // Tube 2 particles (blue) - start near tube2
        particles.push_back({
            tube2_x + spread * uniform_dist(gen),
            tube2_y + spread * uniform_dist(gen),
            0.05 + 0.1 * uniform_dist(gen), // 50-150mm height
            -0.1 + 0.05 * randn(), // Initial velocity toward center
            0.1 + 0.05 * randn(),
            0.02 * randn(),
            1
        });
    }

    std::cout << "Initialized " << particles.size() << " particles\n";

    // Create geometry files for tank and pipes
    createGeometry();

    // Open gnuplot pipe
    FILE* gp = nullptr;
    #ifdef _WIN32
        gp = _popen("gnuplot", "w");
    #else
        gp = popen("gnuplot", "w");
    #endif

    if (gp) {
        // Setup gnuplot for animation with transparency
        fprintf(gp,
            "set terminal qt size 1000,800\n"
            "set xlabel 'X (mm)'\nset ylabel 'Y (mm)'\nset zlabel 'Z (mm)'\n"
            "set xrange [0:400]\nset yrange [0:400]\nset zrange [0:500]\n"
            "set view 60,30\n"
            "set grid\n"
            "set palette defined (0 'red', 1 'blue')\n"
            "set style data points\n"
            "set pointsize 2\n"
            "set ticslevel 0\n"
            "# Set up transparency and line styles\n"
            "set style line 1 lc rgb 'gray' lw 2 dt 1\n"
            "set style line 2 lc rgb 'orange' lw 3 dt 1\n"
            "set style line 3 lc rgb 'green' lw 3 dt 1\n");
        fflush(gp);
        std::cout << "Gnuplot initialized for real-time visualization with geometry\n";
    } else {
        std::cout << "Warning: Could not open gnuplot. Continuing without real-time plots.\n";
    }

    std::cout << "Starting simulation...\n";

    // Simulation loop
    for (int step = 0; step < total_steps; step++) {
        for (auto &p : particles) {
            Vec3 vfield = velocityField(p.x, p.y, p.z);

            double m = (p.source == 0) ? particle_mass1 : particle_mass2;

            // Stokes drag force
            double k_drag = 6 * pi * mu_air * particle_radius / m;
            
            // Apply drag and field velocity
            p.vx += -k_drag * (p.vx - vfield.vx) * dt;
            p.vy += -k_drag * (p.vy - vfield.vy) * dt;
            p.vz += -k_drag * (p.vz - vfield.vz) * dt;

            // Add gravity
            double g = 9.81;
            p.vz -= g * dt;

            // Brownian motion
            double diffusion_factor = sigma_diffusion * sqrt(dt);
            p.vx += diffusion_factor * randn();
            p.vy += diffusion_factor * randn();
            p.vz += diffusion_factor * randn();

            // Update position
            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;

            // Tank boundaries with reflection
            if (p.x <= 0 || p.x >= tank_diam) {
                p.vx = -0.5 * p.vx;
                p.x = std::min(std::max(p.x, 0.001), tank_diam - 0.001);
            }
            if (p.y <= 0 || p.y >= tank_diam) {
                p.vy = -0.5 * p.vy;
                p.y = std::min(std::max(p.y, 0.001), tank_diam - 0.001);
            }
            if (p.z <= 0 || p.z >= tank_height) {
                p.vz = -0.5 * p.vz;
                p.z = std::min(std::max(p.z, 0.001), tank_height - 0.001);
            }
        }

        handleCollisions(particles, particle_radius);

        // Update visualization every 10 steps
        if (step % 10 == 0) {
            // Export to single file that gets overwritten
            exportParticles(particles, "particles.dat");
            
            if (gp) {
                fprintf(gp, "set title '3D Particle Mixing - Step %d/%d'\n", step, total_steps);
                fprintf(gp, 
                    "splot 'tank.dat' with lines ls 1 title 'Tank', \\\n"
                    "      'pipe1.dat' with lines ls 2 title 'Pipe 1 (30°)', \\\n"
                    "      'pipe2.dat' with lines ls 3 title 'Pipe 2 (60°)', \\\n"
                    "      'particles.dat' using 1:2:3:4 with points pt 7 ps 2 lc palette title 'Particles'\n");
                fflush(gp);
            }
            
            if (step % 100 == 0) {
                std::cout << "Step " << step << "/" << total_steps << " completed\n";
            }
        }
        
        // Small delay for smoother animation
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    if (gp) {
        #ifdef _WIN32
            _pclose(gp);
        #else
            pclose(gp);
        #endif
    }

    // Final export
    exportParticles(particles, "particles_final.dat");

    // Compute mixture uniformity
    int bins = 10;
    std::vector<int> bin_count(bins * bins * bins, 0);
    
    for (const auto &p : particles) {
        int ix = std::min(int(p.x / (tank_diam / bins)), bins - 1);
        int iy = std::min(int(p.y / (tank_diam / bins)), bins - 1);
        int iz = std::min(int(p.z / (tank_height / bins)), bins - 1);
        bin_count[ix + iy * bins + iz * bins * bins]++;
    }
    
    double mean = particles.size() / double(bins * bins * bins);
    double variance = 0;
    for (int c : bin_count) {
        variance += (c - mean) * (c - mean);
    }
    variance /= (bins * bins * bins);
    
    double uniformity = (mean > 0) ? 1.0 - sqrt(variance) / mean : 0.0;
    
    std::cout << "\nSimulation completed!\n";
    std::cout << "Final mixture uniformity: " << uniformity << std::endl;
    std::cout << "Final state saved to 'particles_final.dat'\n";
    std::cout << "Geometry files: tank.dat, pipe1.dat, pipe2.dat\n";
    std::cout << "\nTo view final result:\n";
    std::cout << "gnuplot -e \"splot 'tank.dat' with lines lc rgb 'gray', 'pipe1.dat' with lines lc rgb 'orange', 'pipe2.dat' with lines lc rgb 'green', 'particles_final.dat' using 1:2:3:4 with points pt 7 ps 2 lc palette; pause -1\"\n";

    return 0;
}