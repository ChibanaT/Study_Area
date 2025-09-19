#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "../include/math/expr2cmath.hpp"
#include "../include/math/evaluate.hpp"

using namespace std;


int main() {

    cout << "================Newton-Raphson Modified Method================" << endl;
    cout << "Enter the Function: " << endl;
    string expr;
    string expr1;
    cin >> expr1;
    expr = expr1;
    std::string cmath_expr = expr::to_cmath(expr);
    cout << cmath_expr << endl;
    cout << "Enter the Derivative of the Function: " << endl; 
    string expr2;   
    cin >> expr2;
    expr = expr2;
    std::string cmath_der_expr = expr::to_cmath(expr);
    cout << cmath_der_expr << endl;
    cout << "Enter the Second Derivative of the Function: " << endl; 
    string expr3;   
    cin >> expr3;
    expr = expr3;
    std::string cmath_der2_expr = expr::to_cmath(expr);
    cout << cmath_der_expr << endl;

    double x_initial, x_final, tolerance, x, fx_initial, dydx, d2ydx2;
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

    // Data log criation
    std::ofstream data("newton_raphson_modified_data.dat");
    data << expr << endl;
    data << setw(10) << "iteration" << setw(13) << "x" << setw(17) << "f(x)" << setw(18) << "f'(x)" << setw(15) << "f''(x)" << endl;

    // Plot script creation
    ofstream script("plotnewton_raphson_modified.gp");
    script << "set title '" << expr1 << "'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'y'\n";
    script << "plot 'newton_raphson_modified_data.dat' skip 2 using 1:2 with linespoints title 'x', \\\n"
            << "    'newton_raphson_modified_data.dat' skip 2 using 1:3 with linespoints title 'f(x)', \\\n"
            << "    'newton_raphson_modified_data.dat' skip 2 using 1:4 with linespoints title \"f'(x)\",\\\n"
            << "    'newton_raphson_modified_data.dat' skip 2 using 1:4 with linespoints title \"f''(x)\", \n";           
    script.close();

    while (iter < max_iterations) {
        fx_initial = expr::evaluate(cmath_expr,x_initial);        
        dydx = expr::evaluate(cmath_der_expr,x_initial);
        d2ydx2 = expr::evaluate(cmath_der2_expr,x_initial);
        
        x_final = x_initial - (fx_initial/(dydx - 0.5*(fx_initial*d2ydx2/dydx)));

        cout << "Iteration " << iter + 1 << ": x = " << x_final << endl;

        // Saving data
        data << setw(6) << iter << setw(17) << fixed << setprecision(5) << x_initial << setw(17) << fx_initial << setw(17) << dydx << setw(17) << d2ydx2 << endl;

        // Stop criteria
        if (abs(x_final - x_initial) < tolerance && abs(fx_initial) < tolerance) {
            cout << "Converged with x = " << x_final << endl;
            break;            
        }

        x_initial = x_final;
        iter++;
        
    }

    data.close();

    system("gnuplot -persist plotnewton_raphson_modified.gp");

return 0;

}
