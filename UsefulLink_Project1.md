# Some useful material for the Project 1
## Starting FEniCS
Windows (in cmd not powershell):

    docker run -ti -p 127.0.0.1:8000:8000 -v %cd%:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0

MacOS：

    docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0


## Mesh
FEniCS inbuilt mesh

[4. Built-in meshes — FEniCS Project](https://fenicsproject.org/olddocs/dolfin/1.4.0/python/demo/documented/built-in_meshes/python/documentation.html)
##  Calculate the reaction force
[Computing consistent reaction forces — Numerical tours of continuum mechanics using FEniCS master documentation](https://comet-fenics.readthedocs.io/en/latest/demo/tips_and_tricks/computing_reactions.html)
## Plot results on the nodes and quad points with matplotlib
Tripcolor : Plotting the values on the quad points by providing the mesh

[matplotlib.pyplot.tripcolor — Matplotlib 3.7.0 documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tripcolor.html)

Tricontourf: Plotting the values on the node by providing the mesh

[matplotlib.pyplot.tricontourf — Matplotlib 3.7.0 documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tricontourf.html)
