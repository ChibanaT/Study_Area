#include "include/root.hpp"
#include "include/linear_system.hpp"
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

    } else {
        cout << "Invalid choice\n";
    }

    cout << "===========================================================\n";
    return 0;
}
