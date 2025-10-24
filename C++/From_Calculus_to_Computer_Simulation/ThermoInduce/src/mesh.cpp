#include "Mesh.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace TI {

bool generate_annulus_mesh(const MeshParams& params){
    // units: mm in geometry, output is in mm
    double Ro = params.diameter_mm / 2.0;
    double Ri = Ro - params.thickness_mm;
    if(Ri <= 0.0){
        std::cerr << "Invalid geometry: inner radius <= 0\n";
        return false;
    }
    int Nr = std::max(1, params.Nr);
    int Nt = std::max(8, params.Nt);

    std::ofstream fout(params.out_filename);
    if(!fout){
        std::cerr << "Cannot open output file: " << params.out_filename << "\n";
        return false;
    }
    fout << std::fixed << std::setprecision(6);

    // We'll write blocks separated by blank lines so gnuplot can plot lines per ring.
    // For each radial layer (j from 0..Nr), we write the circle at radius r_j sampled by Nt points.
    // r_j will go from Ri to Ro inclusive (Nr+1 rings if we want nodes), but typical structured mesh nodes are Nr+1 radial nodes.
    // We'll write outer and inner rings and intermediate rings.
    for(int j = 0; j <= Nr; ++j){
        double t = static_cast<double>(j) / static_cast<double>(Nr); // 0..1
        double r = Ri + t * (Ro - Ri);
        for(int i = 0; i <= Nt; ++i){
            double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(Nt);
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            fout << x << " " << y << "\n";
        }
        // blank line for gnuplot block separation
        fout << "\n\n";
    }

    // Additionally write radial lines (connect inner to outer at several angles)
    // We'll append blocks that draw radial spokes to visualize connectivity.
    int spokes = std::min(16, Nt);
    for(int s = 0; s < spokes; ++s){
        double theta = 2.0 * M_PI * static_cast<double>(s) / static_cast<double>(spokes);
        for(int j = 0; j <= Nr; ++j){
            double t = static_cast<double>(j) / static_cast<double>(Nr);
            double r = Ri + t * (Ro - Ri);
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            fout << x << " " << y << "\n";
        }
        fout << "\n\n";
    }

    fout.close();
    std::cout << "Wrote mesh to " << params.out_filename << " (outer R="<<Ro<<" mm, inner R="<<Ri<<" mm)\n";
    return true;
}

} // namespace TI
