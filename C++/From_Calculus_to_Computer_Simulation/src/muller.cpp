#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "../include/math/expr2cmath.hpp"
#include "../include/math/evaluate.hpp"

using namespace std;


int main() {

    cout << "================Muller Method================" << endl;
    cout << "Enter the Function: " << endl;
    string expr;
    string expr1;
    cin >> expr1;
    expr = expr1;
    std::string cmath_expr = expr::to_cmath(expr);
    cout << cmath_expr << endl;      

    double x0, x1, x2, x3, x_final, tolerance, f_x0, f_x1, f_x2, f_x3, a, b, c, h1, h2, d1, d2;
    int max_iterations, iter;

    cout << "Enter the x0 value: " << endl;        
    cin >> x0;
    cout << "Enter the x1 value: " << endl;       
    cin >> x1;
    cout << "Enter the x2 value: " << endl;      
    cin >> x2;
    cout << "x0 = " << x0 << "  x1 = " << x1 << "   x2 = " << x2 << endl;        
    cout << "Enter the Tolerance in percentage: " << endl;
    cin >>  tolerance;
    cout << "Enter the Maximum Number of Iterations: " << endl;
    cin >>  max_iterations;

    iter = 0;    
    tolerance /=100.0;  //Convert to percentage
    x3 = 0;
   

    // Data log criation
    std::ofstream data("muller_data.dat");
    data << expr << endl;
    data << setw(10) << "iteration" << setw(13) << "x0" << setw(17) << "x1" << setw(17) << "x2" << setw(17) << "x3" << setw(17) << "f(x0)" << setw(17) << "f(x1)" << setw(17) << "f(x2)" << setw(17) << "f(x3)" << endl;

    // Plot script creation
    ofstream script("plotmuller.gp");
    script << "set title '" << expr1 << "'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'y'\n";
    script << "plot 'muller_data.dat' skip 2 using 1:2 with linespoints title 'x', \\\n"
            << "    'muller_data.dat' skip 2 using 1:3 with linespoints title 'f(x0)', \\\n"
            << "    'muller_data.dat' skip 2 using 1:4 with linespoints title 'f(x1)',\\\n"
            << "    'muller_data.dat' skip 2 using 1:4 with linespoints title 'f(x2)',\\\n"
            << "    'muller_data.dat' skip 2 using 1:4 with linespoints title 'f(x3)', \n";           
    script.close();

    while (iter < max_iterations) {
        f_x0 = expr::evaluate(cmath_expr,x0);
        f_x1 = expr::evaluate(cmath_expr,x1);
        f_x2 = expr::evaluate(cmath_expr,x2);
        f_x3 = expr::evaluate(cmath_expr,x3);        
        h1 = x1 - x0;
        h2 = x2 - x1;
        d1 = (f_x1 - f_x0) / (x1 - x0);
        d2 = (f_x2 - f_x1) / (x2 - x1);

        a = (d2 - d1) / (h2 + h1);
        b = (a * h2) + d2;
        c = f_x2;

        // verify b signal
        int sgn = (b >= 0) ? 1 : -1;
        
        x_final = x2 + ((-2 * c) / (b + sgn * (sqrt(pow(b,2) - 4*a*c))));

        cout << "Iteration " << iter + 1 << ": x = " << x_final << endl;

        // Saving data
        data << setw(6) << iter << setw(17) << fixed << setprecision(5) << x0 << setw(17) << x1 << setw(17) << x2 << setw(17) << x3 << setw(17) << f_x0 << setw(17) << f_x1 << setw(17) << f_x2 << setw(17) << f_x3 << endl;

        // Stop criteria
        double f_x_final = expr::evaluate(cmath_expr, x_final);
        if (abs(x_final - x2) < tolerance && abs(f_x_final) < tolerance) {
            cout << "Converged with x = " << x_final << endl;
            break;            
        }

        x3 = x_final;
        x0 = x1;
        x1 = x2;
        x2 = x3;
        iter++;
        
    }

    data.close();

    system("gnuplot -persist plotmuller.gp");

return 0;

}
