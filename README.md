
# FEniCS_tutorial_2023
This project contains all materials related to FEniCS tutorial in Advanced computational mechanics class (2023 Spring))
## Start FEniCS from Docker
Windows (in cmd not powershell):

    docker run -ti -p 127.0.0.1:8000:8000 -v %cd%:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0

MacOS：

    docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:2016.2.0

## Some methods for solving Stokes equation with FEM

- [Tutorial](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Tutorial/StokesFEM.md)

## Solving Phase Field Equation by FEniCS (03/20/2023)

- [Tutorial](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Tutorial/PhaseFieldDamage_FEniCS.md)
- [Note](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/Note_03202023.jpg?raw=true)
- [Zoom Recording](https://vanderbilt.zoom.us/rec/share/P6Z4lniGENnVe13MfrxiL13Tzgy6ykm3gYDY1WRBz5IJZppgP705C3BRH7vke0P-.mJmepQ_5fFidVP1B)

## Mesh generation for the phase field problem (03/01/2023)
[Tutorial](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Tutorial/MeshforCrack.md) for mshr and Gmsh

[Slides](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Slides/Slides_03012023.pdf) for the class in 03/01/2023
## Project 1 (02/22/2023)
[Useful links for the Project 1](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Tutorial/UsefulLink_Project1.md)

- FEniCS inbuilt mesh generator
- Calculate the reaction force
- Plot results on the nodes and quad points with matplotlib
## Introduction (02/08/2023)
[Slides](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Slides/Mesh%20generation%20and%20coding%20in%20FEniCS.pdf)
