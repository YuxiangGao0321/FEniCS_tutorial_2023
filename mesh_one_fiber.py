# This program generate FEM mesh for a square plate with a fiber
import gmsh
import sys
import numpy as np
from lxml import etree
import os

Pi = np.pi

mesh_size = 0.01
eps = mesh_size*1e-5
L = 1


cx,cy = 0.5,0.5
R = 0.1

gmsh.initialize()
gmsh.model.add("one_fiber")

# Draw geometry by OpenCASCADE (OCC)
gmsh.model.occ.addRectangle(0,0,0,L,L,1)

cir_index = 11
gmsh.model.occ.addCircle(cx, cy, 0,R,cir_index,0,2*Pi)
gmsh.model.occ.addCurveLoop([cir_index],cir_index)
gmsh.model.occ.addPlaneSurface([cir_index],cir_index)

Tool_list = [(2,cir_index)]

gmsh.model.occ.intersect([(2,1)],Tool_list,
	tag = -1,removeObject = False,removeTool = True)
gmsh.model.occ.fragment([(2,1)],Tool_list)

# synchronize OCC model to Gmsh
gmsh.model.occ.synchronize()

surface_list = [i[1] for i in gmsh.model.getEntitiesInBoundingBox(-eps,-eps,-eps,L+eps,L+eps,eps,2)]
#print(surface_list)
gmsh.model.addPhysicalGroup(2, surface_list)

# Define a mesh field if there are several surfaces
# You can define mesh size or add seed with other Gmsh functions 
# (Please refer Gmsh Python API)
gmsh.model.mesh.field.add("Box", 1000)
gmsh.model.mesh.field.setNumber(1000, "VIn", mesh_size)
gmsh.model.mesh.field.setNumber(1000, "VOut", mesh_size)
gmsh.model.mesh.field.setNumber(1000, "XMin", -eps)
gmsh.model.mesh.field.setNumber(1000, "XMax", L+eps)
gmsh.model.mesh.field.setNumber(1000, "YMin", -eps)
gmsh.model.mesh.field.setNumber(1000, "YMax", L+eps)
gmsh.model.mesh.field.setAsBackgroundMesh(1000)

# generate mesh
# You can call different mesh algorithm here.
gmsh.model.mesh.generate(2)

# Save mesh file (.dat)
# dat_file_name = "one_fiber.dat"
# gmsh.write(dat_file_name)

# Uncomment the following code to enable Gmsh popup a window to show your results.
if '-nopopup' not in sys.argv:
	gmsh.fltk.run()

gmsh.finalize()