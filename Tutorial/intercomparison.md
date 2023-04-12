# Numerical Example
The mesh file and the results calculated by FEniCS can be found [here](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/tree/main/Model_intercomparison).

[Here](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Tutorial/FEniCS2Numpy.md) is the tutorial about converting the FEniCS function variable to a Numpy array and loading results from text file.
## 2D uniaxial tension for plane strain
### Problem setup
In this problem, we solve a 2D linear elasticity plane strain problem in a square domain $\Omega := [0,1]^2$ with first Lamé parameters $\lambda$  = 121.5 GPa and Shear modulus $\mu$ = 80.7 GPa. We apply a 0.005 mm y-direction displacement to the top boundary and use rollers on the left and bottom boundaries.

### Analytical solution
Based on the boundary condition, we know that

$$\begin{align}
u_x(x = 0) &= 0\\\\\\
u_y(y = 0) &= 0\\\\\\
u_y(y = 1) &= 0.005
\end{align}$$

As the problem is linear elasticity for plane strain, we can get
$$
\varepsilon_{yy} =  0.005 \\
\varepsilon_{xx} = -\frac{\nu}{1-\nu}\varepsilon_{yy} \\
\varepsilon_{xy} = 0 \\
\sigma_{xx} = 0 \\
\sigma_{xy} = 0 \\
\sigma_{yy} = \frac{E}{(1+\nu)(1-2\nu)}(\nu \varepsilon_{xx} + (1-\nu) \varepsilon_{yy}) \\
u_x = x\varepsilon_{xx} \\
u_y = y\varepsilon_{yy}
$$
where $E = \mu \frac{3\lambda+2\mu}{\lambda + \mu}$ and $\nu = \frac{\lambda}{2(\lambda + \mu)}$.
###  Results evaluation
The solution of this problem solved by FEniCS with the mesh ([Project1/mesh_2d.xml](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Model_intercomparison/Project1/mesh_2d.xml)) can be found in [Project1/result_2d/Data](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/tree/main/Model_intercomparison/Project1/result_2d/Data). The table below shows the L2 error between the FEniCS solution and analytical solution for displacement and stress at each components.
| **Variable** | $u_x$ | $u_y$|$\sigma_{xx}$| $\sigma_{xy}$ | $\sigma_{yy}$|
|--|--|--|--|--|--|--|
| **FEM L2 error** | 2.29e-12 | 1.45e-11 | 1.33e-14 | 4.17e-15 | 1.81e-14 |
The figure of each variable can be found in [Project1/result_2d/Figures](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/tree/main/Model_intercomparison/Project1/result_2d/Figures).

## 3D uniaxial tension
### Problem setup
In this problem, we solve a 3D linear elasticity problem in a cubic domain $\Omega := [0,1]^3$ with first Lamé parameters $\lambda$  = 121.5 GPa and Shear modulus $\mu$ = 80.7 GPa. We apply a 0.005 mm z-direction displacement to the top surface and use rollers on the left, back and bottom boundaries.

### Analytical solution
Based on the boundary condition, we know that
$$
u_x(x = 0) = 0 \\
u_y(y = 0) = 0 \\
u_z(z = 0) = 0 \\
u_z(z = 1) = 0.005
$$
As the problem is linear elasticity, we can get
$$
\varepsilon_{zz} =  0.005 \\
\varepsilon_{xx} = \varepsilon_{yy} =  -\varepsilon_{zz} \frac{\lambda}{2(\mu + \lambda)}\\
\varepsilon_{xy} = \varepsilon_{xz} = \varepsilon_{yz} = 0 \\
\sigma_{xx} = \sigma_{yy} = 0 \\
\sigma_{xy} = \sigma_{xz} = \sigma_{yz} = 0 \\
\sigma_{zz} = \lambda \varepsilon_{xx} + \lambda \varepsilon_{yy}  + (2\mu+\lambda) \varepsilon_{zz} \\
u_x = x\varepsilon_{xx} \\
u_y = y\varepsilon_{yy} \\
u_z = z\varepsilon_{zz}
$$
###  Results evaluation
The solution of this problem solved by FEniCS with the mesh ([Project1/mesh_3d.xml](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Model_intercomparison/Project1/mesh_3d.xml)) can be found in [Project1/result_3d/Data](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/tree/main/Model_intercomparison/Project1/result_3d/Data). The table below shows the L2 error between the FEniCS solution and analytical solution for displacement and stress at each components.
| **Variable** | $u_x$ | $u_y$|$u_z$|$\sigma_{xx}$| $\sigma_{yy}$ | $\sigma_{zz}$|$\sigma_{xy}$| $\sigma_{xz}$ | $\sigma_{yz}$|
|--|--|--|--|--|--|--|--|--|--|--|
| **FEM L2 error** | 3.65e-12 | 3.18e-12 | 1.11e-11 | 3.15e-14 | 2.93e-14 | 4.79e-14 | 4.23e-15 | 7.72e-15 | 7.79e-15|
The pvd file of each variable can be found in [Project1/result_3d/Figures](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/tree/main/Model_intercomparison/Project1/result_3d/Figures).

## 2D uniaxial tension for the phase field fracture
### Problem setup
To verify the staggered numerical implementation in FEniCS, we consider a simple benchmark study single edge cracked specimen subjected to tensile loading. We solve the phase field fracture problem in a square domain $\Omega := [0,1]^2$ with a initial crack from $(0,0.5)$ to $(0.5,0.5)$. We apply a y-direction displacement to the top boundary and fix the bottom boundary, as shown below.

![Schematic diagram of the specimen and the boundary conditions](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/PhaseFieldDomain.jpg?raw=true)

The domain is discretized to finite element mesh ([Project2/mesh_lc=0.015.xml](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/Model_intercomparison/Project2/mesh_lc=0.015.xml)) by Gmsh with the mesh size equal to one-fourth of the length scale ($l_c$ = 0.015). The meshes near the crack path are refined to increase the accuracy and reduce the computational time, as shown below.

![Finite element mesh for $l_c$ = 0.015](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/mesh_lc=0.015.jpg?raw=true)

To compare the load versus displacement diagram with the result in [Sun et al. (2021)](https://www.sciencedirect.com/science/article/pii/S2352431621000626), the same **linear element** and the same material parameters are applied as shown in the table. 
| **Parameter** | **Value** | **Unit** |
|--|--|--|
| $\lambda$ | 121.5 | $\textnormal{kN/mm}^2$|
| $\mu$ | 80.7 | $\textnormal{kN/mm}^2$|
| $G_c$ | 2.7e-3 | $\textnormal{kN/mm}^2$|
| $l_c$ | 0.0015 | $\textnormal{mm}$|
| $\eta$ | 1e-6 | $\textnormal{kN} \cdot \textnormal{s/m}^2$|
And the same displacement control is also applied with a constant increment $\Delta u = 10^{-5}$ when $u<0.005$ and $\Delta u = 10^{-6}$ after that. The strain energy decomposition scheme proposed by Miehe el al. is applied for updating history variable. The iteration will be stopped when the reaction force on the top boundary less than 1% peak load.

### Results
The figure below shows the load versus displacement diagram for the top boundary.
![ the load versus displacement diagram for the top boundary](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/Load-disp%28lc=0.0015%29.jpg?raw=true)
The load displacement data can be found in Project2/result/load_disp_lc=0.015.txt and the y-component of the displacement field and the phase field damage data can be found in Project2/result/Data. The figures of the damage field can be found in Project2/result/Figures.
