#include "Mesh.h"
#include <iostream>

int main(int argc, char** argv){
    TI::MeshParams mp;
    // default geometry: diameter 100 mm, wall 1 mm
    // optional: read simple args (not required). You can edit mp here.
    mp.diameter_mm = 100.0;
    mp.thickness_mm = 1.0;
    mp.Nr = 8;    // radial subdivisions (increase to refine)
    mp.Nt = 180;  // angular samples for smooth circles
    mp.out_filename = "mesh.dat";

    std::cout << "Generating annulus mesh: D=" << mp.diameter_mm << " mm, thickness=" << mp.thickness_mm << " mm\n";
    bool ok = TI::generate_annulus_mesh(mp);
    if(!ok) return 1;

    std::cout << "To plot: gnuplot -persist plots/plot_mesh.plt\n";
    return 0;
}
