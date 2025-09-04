#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "../include/math/expr2cmath.hpp"
#include "../include/math/evaluate.hpp"

using namespace std;


int main() {

    cout << "================Bisection Method================" << endl;
    cout << "Enter the Function: " << endl;
    string expr;
    cin >> expr;
    std::string cmath_expr = expr::to_cmath(expr);
    cout << cmath_expr << endl;

    double x_upper, x_lower, tolerance, x, fx_upper, fx_lower, x_mid, f_mid;
    int max_iterations, iter;

    cout << "Enter the Interval: " << endl;
    cout << "x lower: ";
    cin >>  x_lower;
    cout << "x upper: ";
    cin >>  x_upper;    
    cout << "Enter the Tolerance in percentage: " << endl;
    cin >>  tolerance;
    cout << "Enter the Maximum Number of Iterations: " << endl;
    cin >>  max_iterations;

    iter = 0;    
    tolerance /=100.0;  //Convert to percentage

    
    fx_lower = expr::evaluate(cmath_expr,x_lower);
    fx_upper = expr::evaluate(cmath_expr,x_upper);

    // Data log criation
    std::ofstream data("bisection_data.dat");
    data << expr << endl;
    data << setw(10) << "iteration" << setw(13) << "x lower" << setw(17) << "x upper" << setw(18) << "x middle" << setw(15) << "f(x)" << endl;

    // Plot script creation
    ofstream script("plotbisection.gp");
    script << "set title '" << expr << "'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'y'\n";
    script << "plot 'bisection_data.dat' skip 2 using 1:2 with linespoints title 'x Lower', \\\n"
            << "    'bisection_data.dat' skip 2 using 1:3 with linespoints title 'x Upper', \\\n"
            << "    'bisection_data.dat' skip 2 using 1:4 with linespoints title 'x Middle', \\\n"
            << "    'bisection_data.dat' skip 2 using 1:5 with linespoints title 'f(x)'\n";
    script.close();

    // Signal verification
    if ( fx_lower * fx_upper < 0 ){
 
        while ((iter < max_iterations && abs(x_upper - x_lower) > tolerance)) {  
            x_mid = (x_upper + x_lower) / 2.0;            
            f_mid = expr::evaluate(cmath_expr,x_mid);

            // Saving data
            data << setw(6) << iter << setw(17) << fixed << setprecision(5) << x_lower << setw(17) << x_upper << setw(17) << x_mid << setw(17) << f_mid << endl;
                        
            if (fx_lower * f_mid < 0){
                x_upper = x_mid;
                fx_upper = f_mid;
            }
            
            else if (fx_upper * f_mid < 0){
                x_lower = x_mid;
                fx_lower = f_mid;
            }
            
            cout << "Iteration " << iter + 1 << ": x_mid = " << x_mid << endl;
                        
            iter++;

        }

        cout << "Approximate root after " << iter << " iterations is " << x_mid << endl;  
    }

    
    else {
        cout << " Invalid Interval" << endl;
    }
    data.close();

    system("gnuplot -persist plotbisection.gp");
        
    return 0;

    
}