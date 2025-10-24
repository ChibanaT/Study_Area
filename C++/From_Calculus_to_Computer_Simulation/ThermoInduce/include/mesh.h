#ifndef THERMOINDUCE_MESH_H
#define THERMOINDUCE_MESH_H

#include <string>
#include <vector>

namespace TI {

struct MeshParams {
    double diameter_mm = 100.0;   // outer diameter in mm
    double thickness_mm = 1.0;    // wall thickness in mm
    int Nr = 4;    // radial divisions across wall (small by default)
    int Nt = 64;   // angular divisions
    std::string out_filename = "mesh.dat";
};

// generate annulus structured mesh and write to file suitable for gnuplot
// returns true on success
bool generate_annulus_mesh(const MeshParams& params);

} // namespace TI

#endif // THERMOINDUCE_MESH_H
