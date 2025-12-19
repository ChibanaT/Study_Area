#include "include/root.hpp"
#include "include/linear_system.hpp"
#include "include/interpolation.hpp"
#include "include/ode_solver.hpp"
#include "include/optimization.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace std;
using namespace numerical;

// Function to clear input buffer
void clear_input() {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

int main() {
    cout << "================ Numerical Methods Software ================\n";
    cout << "Choose an option:\n";
    cout << "1. Solve a linear system\n";
    cout << "2. Find roots of a polynomial\n";
    cout << "3. Interpolation\n";
    cout << "4. ODE Solver\n";
    cout << "5. Optimization\n";
    int choice;
    cin >> choice;

    if (choice == 1) {
        // ---------------- Linear System ----------------
        int n;
        cout << "Enter the number of variables: ";
        cin >> n;

        if (n <= 0) {
            cout << "Invalid size\n";
            return 0;
        }

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);

        cout << "Enter the coefficients of the matrix A row by row:\n";
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                cin >> A[i][j];

        cout << "Enter the constants vector b:\n";
        for (int i = 0; i < n; ++i)
            cin >> b[i];

        cout << "Choose a method to solve the linear system:\n";
        cout << "1. Cramer's Rule\n";
        cout << "2. Gaussian Elimination\n";
        cout << "3. Pivoted Gaussian Elimination\n";
        cout << "4. Gauss-Jordan\n";
        cout << "5. LU Decomposition\n";
        cout << "6. Jacobi Iteration\n";
        cout << "7. Gauss-Seidel\n";
        cout << "8. Successive Over-Relaxation (SOR)\n";
        cout << "9. Conjugate Gradient\n";
        cout << "10. Thomas Algorithm (tridiagonal)\n";
        cout << "11. Cholesky Decomposition\n";

        int method;
        cin >> method;

        LinearResult res;

        try {
            switch (method) {
                case 1: res = cramer(A, b); break;
                case 2: res = gaussian_elimination(A, b); break;
                case 3: res = pivot_gauss(A, b); break;
                case 4: res = gauss_jordan(A, b); break;
                case 5: res = lu_decomposition(A, b); break;
                case 6: {
                    int max_iter; double tol;
                    cout << "Enter maximum iterations: "; cin >> max_iter;
                    cout << "Enter tolerance: "; cin >> tol;
                    res = jacobi(A, b, max_iter, tol);
                    break;
                }
                case 7: {
                    int max_iter; double tol;
                    cout << "Enter maximum iterations: "; cin >> max_iter;
                    cout << "Enter tolerance: "; cin >> tol;
                    res = gauss_seidel(A, b, max_iter, tol);
                    break;
                }
                case 8: {
                    int max_iter; double tol; double omega;
                    cout << "Enter maximum iterations: "; cin >> max_iter;
                    cout << "Enter tolerance: "; cin >> tol;
                    cout << "Enter relaxation factor omega (0 < omega < 2): "; cin >> omega;
                    res = sor(A, b, max_iter, tol, omega);
                    break;
                }
                case 9: {
                    int max_iter; double tol;
                    cout << "Enter maximum iterations: "; cin >> max_iter;
                    cout << "Enter tolerance: "; cin >> tol;
                    res = conjugate_gradient(A, b, max_iter, tol);
                    break;
                }
                case 10: res = thomas_algorithm(A, b); break;
                case 11: res = cholesky_decomposition(A, b); break;
                default: throw invalid_argument("Invalid method selection");
            }

            cout << "\n================ Linear System Result ================\n";
            cout << "Converged: " << (res.converged ? "Yes" : "No") << "\n";
            cout << "Iterations: " << res.iterations << "\n";
            cout << "CPU time (s): " << res.cpu_time_sec << "\n";
            cout << "FLOP count: " << res.flop_count << "\n";
            cout << "Solution:\n";
            for (size_t i = 0; i < res.solution.size(); ++i)
                cout << "x[" << i << "] = " << res.solution[i] << "\n";
            cout << "====================================================\n";

        } catch (const exception &e) {
            cout << "Error: " << e.what() << "\n";
        }

    } else if (choice == 2) {
        // ---------------- Polynomial Roots ----------------
        int deg;
        cout << "Enter the degree of the polynomial: ";
        cin >> deg;

        if (deg <= 0) {
            cout << "Degree must be positive\n";
            return 0;
        }

        vector<double> coeffs(deg + 1);
        cout << "Enter the coefficients from highest degree to constant term:\n";
        for (int i = 0; i <= deg; ++i)
            cin >> coeffs[i];

        cout << "Choose a root-finding method:\n";
        cout << "1. Bisection\n2. False Position\n3. Modified False Position\n4. Secant\n";
        cout << "5. Newton-Raphson\n6. Modified Newton\n7. Halley\n8. Brent\n9. Bairstow\n10. Muller Method\n";

        int method;
        cin >> method;

        double tol; int max_iter;
        cout << "Enter tolerance: "; cin >> tol;
        cout << "Enter maximum iterations: "; cin >> max_iter;

        try {
            RootResult res;

            switch (method) {
                case 1: {
                    double a, b; cout << "Enter interval [a b]: "; cin >> a >> b;
                    res = bisect( [&](double x){ return numerical::horner(coeffs, x); }, a, b, tol, max_iter);
                    break;
                }
                case 2: {
                    double a, b; cout << "Enter interval [a b]: "; cin >> a >> b;
                    res = false_position( [&](double x){ return numerical::horner(coeffs, x); }, a, b, tol, max_iter);
                    break;
                }
                case 3: {
                    double a, b; cout << "Enter interval [a b]: "; cin >> a >> b;
                    res = modified_false_position( [&](double x){ return numerical::horner(coeffs, x); }, a, b, tol, max_iter);
                    break;
                }
                case 4: {
                    double a, b; cout << "Enter initial guesses [a b]: "; cin >> a >> b;
                    res = secant( [&](double x){ return numerical::horner(coeffs, x); }, a, b, tol, max_iter);
                    break;
                }
                case 5: {
                    double x0; cout << "Enter initial guess x0: "; cin >> x0;
                    res = newton_raphson( [&](double x){ return numerical::horner(coeffs, x); },
                                         [&](double x){ return numerical::derivative(coeffs, x); },
                                         x0, tol, max_iter);
                    break;
                }
                case 6: {
                    double x0; cout << "Enter initial guess x0: "; cin >> x0;
                    res = modified_newton( [&](double x){ return numerical::horner(coeffs, x); },
                                           [&](double x){ return numerical::derivative(coeffs, x); },
                                           x0, tol, max_iter);
                    break;
                }
                case 7: {
                    double x0; cout << "Enter initial guess x0: "; cin >> x0;
                    res = halley( [&](double x){ return numerical::horner(coeffs, x); },
                                  [&](double x){ return numerical::derivative(coeffs, x); },
                                  [&](double x){ return numerical::second_derivative(coeffs, x); },
                                  x0, tol, max_iter);
                    break;
                }
                case 8: {
                    double a, b; cout << "Enter interval [a b]: "; cin >> a >> b;
                    res = brent( [&](double x){ return numerical::horner(coeffs, x); }, a, b, tol, max_iter);
                    break;
                }
                case 9: {
                    double r, s; cout << "Enter initial guesses for Bairstow (r s): "; cin >> r >> s;
                    res = bairstow(coeffs, r, s, tol, max_iter);
                    break;
                }
                case 10: {
                    double x0, x1, x2;
                    cout << "Enter three initial guesses x0 x1 x2: "; cin >> x0 >> x1 >> x2;
                    res = muller( [&](double x){ return numerical::horner(coeffs, x); }, x0, x1, x2, tol, max_iter);
                    break;
                }
                default: throw invalid_argument("Invalid method selection");
            }

            // Print root-finding results
            cout << "\n================ Root-Finding Result ================\n";
            cout << "Converged: " << (res.converged ? "Yes" : "No") << "\n";
            cout << "Root (last approximation): " << res.root << "\n";
            cout << "Iterations: " << res.iterations << "\n";
            cout << "CPU time (s): " << res.cpu_time_sec << "\n";
            cout << "FLOP count: " << res.flop_count << "\n";
            cout << "History of approximations:\n";
            for (size_t i = 0; i < res.history.size(); ++i)
                cout << "Iter " << i+1 << ": " << res.history[i] << "\n";
            cout << "====================================================\n";

        } catch (const exception &e) {
            cout << "Error: " << e.what() << "\n";
        }

    } else if (choice == 3) {
        // ---------------- Interpolation ----------------
        int n;
        cout << "Enter the number of data points: ";
        cin >> n;

        if (n < 2) {
            cout << "At least 2 points required\n";
            return 0;
        }

        vector<double> x_data(n), y_data(n);
        cout << "Enter x values: ";
        for (int i = 0; i < n; ++i) cin >> x_data[i];
        cout << "Enter y values: ";
        for (int i = 0; i < n; ++i) cin >> y_data[i];

        int n_query;
        cout << "Enter number of query points: ";
        cin >> n_query;
        vector<double> x_query(n_query);
        cout << "Enter query x values: ";
        for (int i = 0; i < n_query; ++i) cin >> x_query[i];

        cout << "Choose interpolation method:\n";
        cout << "1. Linear Interpolation\n";
        cout << "2. Lagrange Polynomial\n";
        cout << "3. Newton Divided Difference\n";
        cout << "4. Cubic Spline\n";

        int method;
        cin >> method;

        try {
            InterpolationResult res;
            double left_deriv = 0.0, right_deriv = 0.0;

            switch (method) {
                case 1: res = linear_interpolation(x_data, y_data, x_query); break;
                case 2: res = lagrange_polynomial(x_data, y_data, x_query); break;
                case 3: res = newton_divided_difference(x_data, y_data, x_query); break;
                case 4: {
                    cout << "Enter left boundary derivative (0 for natural): ";
                    cin >> left_deriv;
                    cout << "Enter right boundary derivative (0 for natural): ";
                    cin >> right_deriv;
                    res = cubic_spline(x_data, y_data, x_query, left_deriv, right_deriv);
                    break;
                }
                default: throw invalid_argument("Invalid method selection");
            }

            cout << "\n================ Interpolation Result ================\n";
            cout << "Success: " << (res.success ? "Yes" : "No") << "\n";
            cout << "CPU time (s): " << res.cpu_time_sec << "\n";
            cout << "FLOP count: " << res.flop_count << "\n";
            cout << "Interpolated values:\n";
            for (size_t i = 0; i < res.interpolated_values.size(); ++i)
                cout << "f(" << x_query[i] << ") = " << res.interpolated_values[i] << "\n";
            cout << "====================================================\n";

        } catch (const exception &e) {
            cout << "Error: " << e.what() << "\n";
        }

    } else if (choice == 4) {
        // ---------------- ODE Solver ----------------
        cout << "Choose ODE solver method:\n";
        cout << "1. Euler\n";
        cout << "2. Runge-Kutta 4th Order (RK4)\n";
        cout << "3. Runge-Kutta-Fehlberg 4(5) (RK45)\n";
        cout << "4. Runge-Kutta 7(8) (RK78)\n";
        cout << "5. Backward Differentiation Formula (BDF)\n";
        cout << "6. Radau\n";

        int method;
        cin >> method;

        double t0, y0, t_end, h;
        cout << "Enter initial time t0: "; cin >> t0;
        cout << "Enter initial value y0: "; cin >> y0;
        cout << "Enter final time t_end: "; cin >> t_end;
        cout << "Enter step size h: "; cin >> h;

        // Simple ODE: dy/dt = -y (exponential decay)
        ODEFunc f = [](double t, double y) { return -y; };

        try {
            ODEResult res;
            double tolerance = 1e-6;

            switch (method) {
                case 1: res = euler(f, t0, y0, t_end, h); break;
                case 2: res = rk4(f, t0, y0, t_end, h); break;
                case 3: res = rk45(f, t0, y0, t_end, h, tolerance); break;
                case 4: res = rk78(f, t0, y0, t_end, h, tolerance); break;
                case 5: res = bdf(f, t0, y0, t_end, h, 2); break;
                case 6: res = radau(f, t0, y0, t_end, h, tolerance); break;
                default: throw invalid_argument("Invalid method selection");
            }

            cout << "\n================ ODE Solver Result ================\n";
            cout << "Converged: " << (res.converged ? "Yes" : "No") << "\n";
            cout << "Steps: " << res.steps << "\n";
            cout << "CPU time (s): " << res.cpu_time_sec << "\n";
            cout << "FLOP count: " << res.flop_count << "\n";
            cout << "Max error: " << res.max_error << "\n";
            cout << "Final value: y(" << t_end << ") = " << res.solution_single.back() << "\n";
            cout << "====================================================\n";

        } catch (const exception &e) {
            cout << "Error: " << e.what() << "\n";
        }

    } else if (choice == 5) {
        // ---------------- Optimization ----------------
        cout << "Choose optimization method:\n";
        cout << "1. Steepest Descent\n";
        cout << "2. BFGS\n";
        cout << "3. Simplex (Linear Programming)\n";
        cout << "4. Sequential Quadratic Programming (SQP)\n";
        cout << "5. Lagrange Multipliers\n";
        cout << "6. Genetic Algorithm\n";
        cout << "7. Particle Swarm Optimization (PSO)\n";

        int method;
        cin >> method;

        int dim;
        cout << "Enter dimension: ";
        cin >> dim;

        vector<double> x0(dim);
        cout << "Enter initial point: ";
        for (int i = 0; i < dim; ++i) cin >> x0[i];

        // Simple objective: f(x) = sum(x_i^2)
        ObjectiveFunc f = [](const vector<double>& x) {
            double sum = 0.0;
            for (double xi : x) sum += xi * xi;
            return sum;
        };

        GradientFunc grad_f = [](const vector<double>& x) {
            vector<double> grad(x.size());
            for (size_t i = 0; i < x.size(); ++i) grad[i] = 2.0 * x[i];
            return grad;
        };

        try {
            OptimizationResult res;
            double tolerance = 1e-6;
            int max_iter = 1000;

            switch (method) {
                case 1: {
                    double step_size;
                    cout << "Enter step size: "; cin >> step_size;
                    res = steepest_descent(f, grad_f, x0, step_size, tolerance, max_iter);
                    break;
                }
                case 2: res = bfgs(f, grad_f, x0, tolerance, max_iter); break;
                case 3: {
                    vector<double> c(dim, 1.0);
                    vector<vector<double>> A;
                    vector<double> b;
                    vector<string> constraint_type;
                    res = simplex(c, A, b, constraint_type, x0);
                    break;
                }
                case 4: {
                    vector<ConstraintFunc> constraints;
                    vector<GradientFunc> grad_constraints;
                    res = sqp(f, grad_f, constraints, grad_constraints, x0, tolerance, max_iter);
                    break;
                }
                case 5: {
                    vector<ConstraintFunc> equality_constraints;
                    vector<GradientFunc> grad_equality_constraints;
                    res = lagrange_multipliers(f, grad_f, equality_constraints, grad_equality_constraints, x0, tolerance, max_iter);
                    break;
                }
                case 6: {
                    vector<double> lower_bounds(dim, -10.0);
                    vector<double> upper_bounds(dim, 10.0);
                    res = genetic_algorithm(f, dim, lower_bounds, upper_bounds, 50, 100, 0.1, 0.8);
                    break;
                }
                case 7: {
                    vector<double> lower_bounds(dim, -10.0);
                    vector<double> upper_bounds(dim, 10.0);
                    res = pso(f, dim, lower_bounds, upper_bounds, 30, 100, 0.7, 1.5, 1.5);
                    break;
                }
                default: throw invalid_argument("Invalid method selection");
            }

            cout << "\n================ Optimization Result ================\n";
            cout << "Converged: " << (res.converged ? "Yes" : "No") << "\n";
            cout << "Iterations: " << res.iterations << "\n";
            cout << "CPU time (s): " << res.cpu_time_sec << "\n";
            cout << "FLOP count: " << res.flop_count << "\n";
            cout << "Optimal value: " << res.optimal_value << "\n";
            cout << "Optimal point: [";
            for (size_t i = 0; i < res.optimal_point.size(); ++i) {
                cout << res.optimal_point[i];
                if (i < res.optimal_point.size() - 1) cout << ", ";
            }
            cout << "]\n";
            cout << "====================================================\n";

        } catch (const exception &e) {
            cout << "Error: " << e.what() << "\n";
        }

    } else {
        cout << "Invalid choice\n";
    }

    cout << "===========================================================\n";
    return 0;
}
