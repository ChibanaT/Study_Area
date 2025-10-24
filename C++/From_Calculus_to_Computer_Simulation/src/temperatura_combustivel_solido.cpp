# include <iostream>
# include <vector>
# include <cmath>

using namespace std;

    //////////////////////////////////////////////////////////////////////
    //                         Exact Solution                          //
   /////////////////////////////////////////////////////////////////////


//Transcendental function
double f(double lambda, double Bi) {
    return tan(lambda) - lambda / Bi;
}

//Bisection method
double bisection(double a, double b, double Bi, double tol = 1e-6) {
    double fa = f(a, Bi);
    double fb = f(b, Bi);    
    if (f(a, Bi) * f(b, Bi) > 0) {
        cerr << "Invalid interval: function does not change sign." << endl;
        return -1;
    }
    double c;
    while ((b - a) / 2 > tol){
        c = (a + b) / 2;
        double fc = f(c, Bi);
        if (fabs(fc) < tol) return c;
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return (a + b) / 2;
}

//Analytical steady-state solution
vector<double> exact_solution(int n, double L, double Bi, double T_inf, double T_0) {
    vector<double> T(n);
    double lambda = bisection(0.01, 1.4, Bi);       //Avoid singularity

    double denom = cos(lambda) - (sin(lambda) / lambda) * Bi;

    for (int i = 0; i < n; ++i) {
        double x = L * i/ (n - 1);
        double numerator = cos(lambda * x / L) - (sin(lambda * x / L) / lambda) * Bi;
        T[i] = (numerator / denom) * (T_inf - T_0) + T_0;
    }
    return T;
}

    //////////////////////////////////////////////////////////////////////
    //                        Thomas Algorithm                         //
   /////////////////////////////////////////////////////////////////////

vector<double> thomas_algorithm(
    const vector<double>& a,
    const vector<double>& b,
    const vector<double>& c,
    const vector<double>& d
) {
    int n = b.size();
    vector<double> c_star(n, 0.0);
    vector<double> d_star(n, 0.0);
    vector<double> x(n);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for(int i = 1; i < n; ++i){
        double m = b[i] - a[i - 1] * c_star[i - 1];
        if (i < n - 1)
            c_star[i] = c[i] / m;
        d_star[i] = (d[i] - a[i - 1] * d_star[i - 1]) / m;
    }

    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    return x;
}

int main(){

    //////////////////////////////////////////////////////////////////////
    //                         Creating Matrix                         //
   /////////////////////////////////////////////////////////////////////

    //Matrix dimension
    int n = 7;

    //Variables
    double A = 0.0;
    double T_0 = 350.0;
    double T_inf = 300.0;
    double Q = 10.0;
    double h = 25.0;
    double dx = 0.1;
    double k = 0.5;
    double F0 = 0.1;
    double Bi = h * dx / k;

    //Diagonal vectors:
    vector<double> a(n-1);      //Subdiagonal
    vector<double> b(n);        //Diagonal
    vector<double> c(n-1);       //Superdiagonal

    //Fill diagonal (a)
    a[0] = -F0;
    a[1] = -F0;
    a[2] = -F0;
    a[3] = -F0;
    a[4] = -F0;
    a[5] = -2 * F0;

    //Fill principal diagonal (b)
    b[0] = 1 + 2 * F0;
    b[1] = 1 + 2 * F0;
    b[2] = 1 + 2 * F0;
    b[3] = 1 + 2 * F0;
    b[4] = 1 + 2 * F0;
    b[5] = 1 + 2 * F0;
    b[6] = 1 + 2 * F0 + 2 * F0 * Bi ;

    //Fill superdiagonal (c)
    c[0] = -2 * F0;
    c[1] = -F0;
    c[2] = -F0;
    c[3] = -F0;
    c[4] = -F0;
    c[5] = -F0;

    // temperature
    vector<double> T_p{300, 300, 300, 300, 300, 300, 300};

    //Independent term
    vector<double> B(n, 0.0);
    for (int i = 0; i < n-1; ++i)
        B[i] = T_p[i] + A;

    B[n-1] = T_p[n-1] + F0 * (2 * Bi * T_inf) + Q * dx * dx / k;

    //Using Thomas algorithm
    vector<double> T_next = thomas_algorithm(a, b, c, B);

    //Print Results
    cout << "Future Temperatures (Thomas Method): " << endl;
    for (double temp : T_next) {
        cout << temp << "   ";
    }
    cout << endl;

    //Exact solution (A = 0)
    if (A == 0.0) {
        vector<double> T_exact = exact_solution(n, dx * (n - 1), Bi, T_inf, T_p[0]);
        cout << "Exact Solution (A = 0): " << endl;
        for (double temp : T_exact) {
            cout << temp << "    ";
        }        
    }
    cout << endl;

    return 0;
}
    