#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

/*
 Coupled radiative–convective heater model
 Convection coefficient computed from Re, Pr, Nu
 References:
 - Incropera et al., Fundamentals of Heat and Mass Transfer
 - Cengel, Heat Transfer: A Practical Approach
 - Chapra, Numerical Methods for Engineers
*/

int main() {

    // -----------------------------
    // Physical constants
    // -----------------------------
    const double sigma = 5.670374419e-8; // Stefan–Boltzmann [W/m2.K4]

    // -----------------------------
    // Air properties at ~200 °C
    // (Incropera – property tables)
    // -----------------------------
    const double rho = 0.746;        // Density [kg/m3]
    const double cp  = 1005.0;       // Specific heat [J/kg.K]
    const double k   = 0.036;        // Thermal conductivity [W/m.K]
    const double mu  = 2.57e-5;      // Dynamic viscosity [Pa.s]

    // -----------------------------
    // System parameters
    // -----------------------------
    const double Pel = 2540.0;       // Electrical power [W]
    const double epsilon = 0.8;      // Emissivity [-]

    const double Ar = 0.0849;        // Resistance area [m2]
    const double At = 0.157;         // Tube inner area [m2]

    const double Dh = 0.096;         // Hydraulic diameter [m]
    const double Dflow = 0.096;      // Flow diameter [m]

    const double L = 0.5;            // Length [m]

    const double Vdot = 2.4 / 60.0;  // Volumetric flow rate [m3/s]

    const double Tin = 25.0 + 273.15;
    const double Tout = 200.0 + 273.15;
    const double Ta = Tout;          // Bulk air temperature [K]

    // -----------------------------
    // Flow calculations
    // -----------------------------
    double Aflow = M_PI * pow(Dflow/2.0, 2);
    double u = Vdot / Aflow;                    // Mean velocity [m/s]

    double Re = rho * u * Dh / mu;
    double Pr = cp * mu / k;

    // Dittus–Boelter (turbulent internal flow)
    double Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4);

    double h = Nu * k / Dh;                     // Convective coefficient

    // -----------------------------
    // Newton–Raphson parameters
    // -----------------------------
    const int maxIter = 100;
    const double tol = 1e-6;
    const double relax = 0.3;

    // Initial guesses (K)
    double Tr = 700.0;
    double Tt = 600.0;

    // -----------------------------
    // Newton–Raphson loop
    // -----------------------------
    for (int iter = 0; iter < maxIter; ++iter) {

        double f1 = Pel
            - epsilon * sigma * Ar * (pow(Tr,4) - pow(Tt,4))
            - h * Ar * (Tr - Ta);

        double f2 = epsilon * sigma * Ar * (pow(Tr,4) - pow(Tt,4))
            - h * At * (Tt - Ta);

        double df1_dTr = -epsilon * sigma * Ar * (4.0 * pow(Tr,3)) - h * Ar;
        double df1_dTt =  epsilon * sigma * Ar * (4.0 * pow(Tt,3));

        double df2_dTr =  epsilon * sigma * Ar * (4.0 * pow(Tr,3));
        double df2_dTt = -epsilon * sigma * Ar * (4.0 * pow(Tt,3)) - h * At;

        double detJ = df1_dTr * df2_dTt - df1_dTt * df2_dTr;
        if (fabs(detJ) < 1e-12) {
            cout << "Jacobian singular\n";
            return -1;
        }

        double dTr = (-f1 * df2_dTt + f2 * df1_dTt) / detJ;
        double dTt = (-df1_dTr * f2 + df2_dTr * f1) / detJ;

        Tr += relax * dTr;
        Tt += relax * dTt;

        if (fabs(dTr) < tol && fabs(dTt) < tol) break;
    }

    // -----------------------------
    // Results
    // -----------------------------
    cout << fixed << setprecision(2);
    cout << "\n===== FLOW DATA =====\n";
    cout << "Velocity: " << u << " m/s\n";
    cout << "Reynolds: " << Re << "\n";
    cout << "Prandtl: " << Pr << "\n";
    cout << "Nusselt: " << Nu << "\n";
    cout << "h: " << h << " W/m2.K\n";

    cout << "\n===== TEMPERATURES =====\n";
    cout << "Resistance: " << Tr - 273.15 << " °C\n";
    cout << "Tube:       " << Tt - 273.15 << " °C\n";
    cout << "Air:        " << Ta - 273.15 << " °C\n";

    return 0;
}
