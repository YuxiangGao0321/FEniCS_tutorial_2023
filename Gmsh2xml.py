import sys
import numpy as np
from lxml import etree
import os

# input the name of the mesh file (.dat)
dat_file_name = "mesh.dat"

node_data = np.zeros((int(1e8),3))
ele_data = np.zeros((int(1e8),4))
ni = 0
ei = 0
with open(dat_file_name) as f:
    for line in f:
        dataline = line.split(" ")
        if dataline[0] == "node":
            node_data[ni,0] = ni
            node_data[ni,1] = float(dataline[2])
            node_data[ni,2] = float(dataline[3])
            ni = ni +1
        elif dataline[0] == "element":

            ele_data[ei,0] = ei
            ele_data[ei,1] = dataline[4]
            ele_data[ei,2] = dataline[5]
            ele_data[ei,3] = dataline[6]
            ei = ei+1
        else:
            continue

nodes = [
    etree.Element('vertex', index=str(int(node_i[0])), x="%.5f" % node_i[1], y ="%.5f" % node_i[2])
    for node_i in node_data[:ni,:]
    ]
eles = [
    etree.Element('tetrahedron',index=str(int(ele_i[0])),v0=str(int(ele_i[1])-1),v1 =str(int(ele_i[2])-1),v2 =str(int(ele_i[3])-1))
    for ele_i in ele_data[:ei,:]
    ]

#===================================================
if not os.path.exists("input"):
    os.mkdir("input")


ns = "http://fenicsproject.org"
dolfin = etree.Element("dolfin", nsmap = {"dolfin":ns})
mesh = etree.Element("mesh",celltype="triangle",dim="2" )
vertices = etree.Element("vertices",size=str(ni))
dolfin.append(mesh)
mesh.append(vertices)
vertices.extend(nodes)
cells = etree.Element("cells",size=str(ei))
mesh.append(cells)
cells.extend(eles)
et = etree.ElementTree(dolfin)
et.write('input/hole_mesh.xml',pretty_print=True,xml_declaration=True)
# np.savetxt("nodes.txt.gz",node_data[:ni,:])
# np.savetxt("elems.txt.gz",ele_data[:ei,:])
print("Finished")

# The following code will remove .dat file
# if os.path.exists(dat_file_name):
# 		os.remove(dat_file_name)
# else:
# 	print(dat_file_name)
# 	print("The file does not exist")