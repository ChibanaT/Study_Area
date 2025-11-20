#include "../../include/root.hpp"
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace numerical {

RootResult brent(const Func &f, double a, double b, double tol, int max_iter) {

    auto start_time = std::chrono::high_resolution_clock::now();

    RootResult result;
    result.iterations = 0;
    result.converged = false;
    result.flop_count = 0;

    double fa = f(a), fb = f(b);
    result.flop_count += 2;

    if (fa == 0.0) {
        result.root = a;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (fb == 0.0) {
        result.root = b;
        result.converged = true;
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        result.cpu_time_sec = elapsed.count();
        return result;
    }
    if (fa * fb > 0.0)
        throw std::invalid_argument("f(a) and f(b) must have opposite signs");

    double c = a, fc = fa;
    double d = c;
    double s = b;
    double fs = fb;
    bool mflag = true; // was last step a bisection?

    for (int k = 0; k < max_iter; ++k) {
        // Ensure |f(a)| < |f(b)|
        if (std::fabs(fc) < std::fabs(fb)) {
            std::swap(a, b); std::swap(fa, fb);
            std::swap(c, a); std::swap(fc, fa);
        }

        // Convergence check
        if (std::fabs(fb) < tol || std::fabs(b - a) < tol) {
            result.root = b;
            result.converged = true;
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            result.cpu_time_sec = elapsed.count();
            return result;
        }

        // Attempt inverse quadratic interpolation or secant
        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = (a*fb*fc)/((fa-fb)*(fa-fc))
              + (b*fa*fc)/((fb-fa)*(fb-fc))
              + (c*fa*fb)/((fc-fa)*(fc-fb));
            result.flop_count += 10; // rough estimate
        } else {
            // Secant method
            if (std::fabs(fb - fa) < 1e-14)
                throw std::runtime_error("Division by zero in Brent's secant step");
            s = b - fb * (b - a) / (fb - fa);
            result.flop_count += 6;
        }

        // Conditions to perform bisection instead
        if (!((s > (3*a + b)/4 && s < b) || (s < (3*a + b)/4 && s > b)) ||
            (mflag && std::fabs(s - b) >= std::fabs(b - c)/2) ||
            (!mflag && std::fabs(s - b) >= std::fabs(c - d)/2)) {
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }

        fs = f(s);
        result.history.push_back(s);
        ++result.iterations;
        result.flop_count += 1;

        // Shift variables
        d = c; c = b; fc = fb;
        if (fa * fs < 0.0) { b = s; fb = fs; } else { a = s; fa = fs; }
        if (std::fabs(fa) < std::fabs(fb)) { std::swap(a, b); std::swap(fa, fb); }
    }

    // If not converged
    result.root = b;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    result.cpu_time_sec = elapsed.count();

    return result;
}

} // namespace numerical
    