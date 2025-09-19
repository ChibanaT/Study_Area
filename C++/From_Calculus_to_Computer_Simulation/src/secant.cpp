# include <iostream>
# include <fstream>
# include <cstdlib>
# include <cmath>
# include <iomanip>
# include "../include/math/expr2cmath.hpp"
# include "../include/math/evaluate.hpp"

using namespace std;

int main(){

    cout << "================Regula Falsi Method================" << endl;
    cout << "Enter the Function: " << endl;
    string expr;
    cin >> expr;
    std::string cmath_expr = expr::to_cmath(expr);
    cout << cmath_expr << endl;

    double x, x_final, x_initial, tolerance, fx_final, fx_initial;
    int max_iterations, iter;

    cout << "Enter the Interval: " << endl;
    cout << "x initial: ";
    cin >> x_initial;
    cout << "x final: ";
    cin >> x_final;
    cout << "Enter the tolerance in percentage: " << endl;
    cin >> tolerance;
    cout << "Enter the Maximum Number of Iterations: " << endl;
    cin >> max_iterations;

    iter = 0;
    tolerance /=100.0; //convert to percentage

    // Data log creation
    std::ofstream data("secant_data.dat");
    data << expr << endl;
    data << setw(10) << "Iteration" << setw(13) << "x initial" << setw(17) << "x final" << setw(15) << "f(x)" << endl;

    // Plot script creation
    ofstream script("plotsecant.gp");
    script << "set title '" << expr << "'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'y'\n";
    script << "plot 'secant_data.dat' skip 2 using 1:2 with linespoints title 'x initial', \\\n"
            << "    'secant_data.dat' skip 2 using 1:3 with linespoints title 'x final', \\\n"
            << "    'secant_data.dat' skip 2 using 1:4 with linespoints title 'f(x)'\n";
    script.close();

    while (iter < max_iterations){
        fx_initial = expr::evaluate(cmath_expr,x_initial);
        fx_final = expr::evaluate(cmath_expr,x_final);

        x = x_initial - fx_initial * ((x_initial - x_final) / (fx_initial - fx_final));

        cout << "Iteration " << iter + 1 << ": x = " << x << endl; 

        // Saving data
        data << setw(6) << iter << setw(17) << fixed << setprecision(5) << x_initial << setw(17) << x_final << setw(17) << fx_final << endl;

        if (abs(fx_final) < tolerance){
            cout << "Converged with x = " << x_final << endl;
            break;
        }

        x_initial = x_final;
        x_final = x;
        iter++;
 
    }

    data.close();+

    system("gnuplot -persist plotsecant.gp");
            
    return 0;
}