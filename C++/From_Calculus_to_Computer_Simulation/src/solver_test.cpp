#include <iostream>
#include <vector>
#include <chrono>
#include "rk_solvers.hpp"     // <- usando seu arquivo

// --------------------- System: 3 reactors (Chapra) ---------------------
// C1' = (C3 - C1)/tau1
// C2' = (C1 - C2)/tau2
// C3' = (C2 - C3)/tau3

struct Params { double tau1,tau2,tau3; };

std::vector<double> reactors(double t, std::vector<double> C, Params p)
{
    return {
        (C[2]-C[0])/p.tau1,
        (C[0]-C[1])/p.tau2,
        (C[1]-C[2])/p.tau3
    };
}

// ----------------------- Solve and compare -----------------------------
int main()
{
    Params pr = {3,2,5};         // time constants of each reactor
    std::vector<double> y0 = {1,0,0};     // initial concentration

    double t0=0, tf=50;

    // -------------------- RK45 --------------------
    RK45 rk45(0.05,1e-6);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto result45 = rk45.solve(
        [&](double t,std::vector<double> y){return reactors(t,y,pr);},
        t0,y0,tf
    );
    auto t2 = std::chrono::high_resolution_clock::now();
    double time45 = std::chrono::duration<double>(t2-t1).count();


    // -------------------- RK78 --------------------
    RK78 rk78(0.05,1e-10);
    auto t3 = std::chrono::high_resolution_clock::now();
    auto result78 = rk78.solve(
        [&](double t,std::vector<double> y){return reactors(t,y,pr);},
        t0,y0,tf
    );
    auto t4 = std::chrono::high_resolution_clock::now();
    double time78 = std::chrono::duration<double>(t4-t3).count();


    // -------------------- Output --------------------
    std::cout << "\nFinal concentrations at t="<<tf<<"\n";
    std::cout << "\nRK45 -> C1="<<result45[0]<<"  C2="<<result45[1]<<"  C3="<<result45[2];
    std::cout << "\nTime RK45 = " << time45 << " s";

    std::cout << "\n\nRK78 -> C1="<<result78[0]<<"  C2="<<result78[1]<<"  C3="<<result78[2];
    std::cout << "\nTime RK78 = " << time78 << " s\n";

    std::cout << "\nDifference RK45 vs RK78:\n";
    std::cout << "dC1="<<fabs(result45[0]-result78[0])
              << "  dC2="<<fabs(result45[1]-result78[1])
              << "  dC3="<<fabs(result45[2]-result78[2])<<"\n\n";
}
