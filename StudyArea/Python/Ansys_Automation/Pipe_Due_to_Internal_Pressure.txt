#  Launch PyMAPDL

import numpy as np
import matplotlib.pyplot as plt
from ansys.mapdl.core import launch_mapdl

# start mapdl
mapdl = launch_mapdl()
print(mapdl)

# Setup the pipe cross section

def pipe_plane_strain(e, nu, inn_radius, out_radius, press, AESIZE):
    
    mapdl.clear()
    mapdl.prep7()
   
    mapdl.et(1, "PLANE182", kop3=2)

    mapdl.pcirc(inn_radius, out_radius, theta1=0, theta2=90)
    mapdl.cm("PIPE_PROFILE", "AREA")

    mapdl.mp("EX", 1, e)
    mapdl.mp("PRXY", 1, nu)

    
    mapdl.aesize("ALL", AESIZE)
    mapdl.mshape(0, "2D") 
    mapdl.mshkey(1) 
    mapdl.cmsel("S", "PIPE_PROFILE")
    mapdl.amesh("ALL")

    mapdl.nsel("S", "LOC", "X", 0)  
    mapdl.cm("X_FIXED", "NODES")

    mapdl.nsel("S", "LOC", "Y", 0)
    mapdl.cm("Y_FIXED", "NODES")
    mapdl.allsel()

    mapdl.lsel("S", "RADIUS", vmin=rad1)
    mapdl.cm("PRESSURE_EDGE", "LINE")
    mapdl.allsel()

    mapdl.slashsolu()
    mapdl.antype("STATIC", "NEW")

    mapdl.d("X_FIXED", "UX", 0)
    mapdl.d("Y_FIXED", "UY", 0)
    
    mapdl.csys(1)
    
    mapdl.sfl("PRESSURE_EDGE", "PRES", press)

    
    mapdl.allsel()
    mapdl.solve()
    mapdl.finish()
    
    mapdl.post1()
    mapdl.set(1, 1)

    max_eqv_stress = np.max(mapdl.post_processing.nodal_eqv_stress())
    all_dof = mapdl.mesh.nnum_all
    num_dof = 2*all_dof.size

    return num_dof, max_eqv_stress

# Perform Mesh Convergence Study

rad1 = 175
rad2 = 200
pressure = 100

e = 2e5
nu = 0.3

num_dof = []
max_stress = []

esizes = np.logspace(1.4, 0, 20)

for esize in esizes:
    dof, eqv_stress = pipe_plane_strain(e, nu, rad1, rad2, pressure, esize)
    num_dof.append(dof)
    max_stress.append(eqv_stress)
    print(f'DOF: {dof:5d}   #Stress: {eqv_stress:.2f} MPa')

# Plot mesh convergence results

plt.figure(figsize=(8, 6), dpi=80)
plt.plot(num_dof, max_stress, 'b-o')
plt.plot([num_dof[0], num_dof[-1]], [max_stress[-1], max_stress[-1]], 'r:')
plt.title('Mesh Convergence Study')
plt.xlabel('Number of DOF')
plt.ylabel('Maximum eqv. Stress (MPa)')
plt.show()

# Plot results from converged mesh analysis

mapdl.allsel('ALL')
mapdl.eplot(
    title='Element Plot', line_width=1, show_bounds=True, cpos="xy"
)

# Plot nodal displacement

mapdl.post1()
mapdl.set(1, 1)

mapdl.post_processing.plot_nodal_displacement(
    'NORM', cpos="xy", cmap="magma",
)

# Plot nodal equivalent stress

mapdl.post_processing.plot_nodal_eqv_stress(cpos="xy", cmap="magma")

# Exit MAPDL

mapdl.exit()