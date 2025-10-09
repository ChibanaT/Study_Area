#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>


namespace numerical {

RootResult bisect(const Func &f, double a, double b, double tol, int max_iter){

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double fa = f(a), fb = f(b);
    result.flop_count += 2;
    if (fa == 0.0){
        result.root = a;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;

    }
    if (fb == 0.0){
        result.root = b;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;

    }
    if (fa * fb > 0.0)
        throw std::invalid_argument("f(a) and f(b) must have opposite signs");

    double c = a;
    for (int k = 0; k < max_iter; ++k){
        c = 0.5 * (a + b);
        result.history.push_back(c);
        double fc = f(c);
        result.flop_count +=3; 
        ++result.iterations;
        if (std::fabs(fc) < tol || 0.5 * (b - a) < tol){
            result.root = c;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }
        if (fa * fc < 0.0){b = c; fb = fc;} else {a = c; fa = fc;};
        result.flop_count +=1;

    }
    result.root = c;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} //namespace numerical