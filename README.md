
# FEniCS_tutorial_2023
This project contains all materials related to FEniCS tutorial in Advanced computational mechanics class (2023 Spring))
## Start FEniCS from Docker
Windows (in cmd not powershell):

    docker run -ti -p 127.0.0.1:8000:8000 -v %cd%:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0

MacOS：

    docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0

## Mesh generation for the phase field problem (03/01/2023)
[Tutorial](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/MeshforCrack.md) for mshr and Gmsh
## Project 1 (02/22/2023)
[Useful links for the Project 1](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/UsefulLink_Project1.md)
## Introduction (02/08/2023)
[Slides](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Mesh%20generation%20and%20coding%20in%20FEniCS.pdf)
