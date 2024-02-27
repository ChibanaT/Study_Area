# All units in (m, Kg, s)

LENGTH = 7
WIDTH = 3
DEPTH = 0.2
RADIUS = 0.5
NUM = 3

E = 2e11
NU = 0.27

PRESSURE = 1000

# launch MAPDL and create geometry

from ansys.mapdl.core import launch_mapdl
mapdl = launch_mapdl()

mapdl.clear()
mapdl.prep7()
mapdl.block(0, LENGTH, 0, WIDTH, 0, DEPTH)
for i in range(1,NUM+1):
    mapdl.cyl4(i*LENGTH/(NUM+1),WIDTH/2,RADIUS,'','','',2*DEPTH)
mapdl.vsbv(1,'all')
mapdl.vplot('all')

# define material properties, mesh attributes and generate mesh

mapdl.lesize("ALL", 0.15, layer1=1)

mapdl.mp('ex',1,E)
mapdl.mp('nuxy',1,NU)

mapdl.et(1,'SOLID186')
mapdl.mshape(1, "3D")
mapdl.mshkey(0)
mapdl.vmesh('all')
mapdl.eplot()

# apply loads and boundary conditions

mapdl.nsel('s','loc','x',0)
mapdl.d('all','all',0)

mapdl.nsel('s','loc','x', LENGTH)
mapdl.sf('all','pres',PRESSURE)

mapdl.allsel()
mapdl.finish()

# solve the static problem

mapdl.slashsolu()
mapdl.solve()
mapdl.finish()

# enter the solver routine and solve 
mapdl.slashsolu()
output = mapdl.solve()

print(output)

# plot the stress contours for the model

result = mapdl.result
result.plot_principal_nodal_stress(0, 'seqv', background = 'w', show_edges = True, text_color = 'k', add_text = True)

# exit MAPDL

mapdl.exit()