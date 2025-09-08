#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "../include/math/expr2cmath.hpp"
#include "../include/math/evaluate.hpp"

using namespace std;


int main() {

    cout << "================Newton-Raphson Method================" << endl;
    cout << "Enter the Function: " << endl;
    string expr;
    cin >> expr;
    std::string cmath_expr = expr::to_cmath(expr);
    cout << cmath_expr << endl;
    cout << "Enter the Derivative of the Function: " << endl;    
    cin >> expr;
    std::string cmath_der_expr = expr::to_cmath(expr);
    cout << cmath_der_expr << endl;

    double x_initial, x_final, tolerance, x, fx_initial, dydx;
    int max_iterations, iter;

    cout << "Enter the Initial Value: " << endl;
    cout << "x initial: ";
    cin >>  x_initial;        
    cout << "Enter the Tolerance in percentage: " << endl;
    cin >>  tolerance;
    cout << "Enter the Maximum Number of Iterations: " << endl;
    cin >>  max_iterations;

    iter = 0;    
    tolerance /=100.0;  //Convert to percentage
    x_final = 0;

    fx_initial = expr::evaluate(cmath_expr,x_initial);
    dydx = expr::evaluate(cmath_der_expr,x_initial);

    while (iter < max_iterations) {
        fx_initial = expr::evaluate(cmath_expr,x_initial);        
        dydx = expr::evaluate(cmath_der_expr,x_initial);
        
        x_final = x_initial - (fx_initial/dydx);

        cout << "Iteration " << iter + 1 << ": x = " << x_final << endl;

        // Stop criteria
        if (abs(x_final - x_initial) < tolerance && abs(fx_initial) < tolerance) {
            cout << "Converged with x = " << x_final << endl;
            break;            
        }

        x_initial = x_final;
        iter++;
        
    }

return 0;

}
