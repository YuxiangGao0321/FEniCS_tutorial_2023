# Solving Phase Field Equation by FEniCS
## Phase Field Equation for Damage Evolution
The strong form of the phase fields equation for damage evolution is

$$\eta \dot{D} = l_c \Delta D - \frac{D}{l_c} + 2(1-D)\frac{\mathcal{H}}{\mathcal{G}_c},$$

where the history variable $\mathcal{H}$ is defined by the max value in the loading history of the positive component of elastic strain energy.

$$\mathcal{H} = {max}_{\tau \in [0,t]}\ [\psi_e^+ (\epsilon( \mathbf{x},\tau))]$$

## Staggered scheme for solving the Phase Field Equation
[Sun et al.](https://www.sciencedirect.com/science/article/abs/pii/S2352431621000626) employed a three-step staggered numerical solution strategy at any given pseudo time increment:

 1. Compute the displacement $u$ using static equilibrium along with its boundary conditions and the constitutive damage law;
 2. Evaluate the tensile elastic strain energy $\psi_e^+$ and/or the history variable $\mathcal{H}$, and then update the crack driving force $\mathcal{H}/\mathcal{G_c}$;
 3. Compute the isotropic damage variable $D$ by solving the phase field evolution equation.

## Defining the Function Space and Variables for the Phase Field Equation
Previously, we discussed solving elasticity equations and calculating the reaction forces in FEniCS.

In FEniCS implementation, similar to the elasticity problem, for solving the phase field equation, a scalar function space on the node for the damage field need to defined by

    p = 1 # the degree of the elements
    S = FunctionSpace(mesh,"Lagrange",p) # phase-field shape function

And a scalar function space on the quad points for the damage field can be defined by

    deg_quad = 1 # the degree of the elements
    scalar_quad = FiniteElement("Quadrature", cell=mesh.ufl_cell(),
	    degree=deg_quad, quad_scheme="default")
    SQ = FunctionSpace(mesh, scalar_quad)  # quadrature points in scalar space
    form_params = {"quadrature_degree": deg_quad}
 
  The damage field and degradation function need to be defined on the quad points and initialized by

    deg = Function(SQ) # degradation function 
    deg.vector()[:] = 1
    dmg = Function(SQ) # damge field
    dmg.vector()[:] = 0

The test and trial functions for the phase field equation need to be defined on the nodes

    d = TrialFunction(S)
    omega = TestFunction(S)

Finally, we define two variables, a Numpy array and a FEniCS function variable, to save the history variable after each iteration.

    kappa = np.zeros(nQ_local)
    phi = Function(SQ)

Where `nQ_local` is the number of the quad points, which can be derived from

    ### Coordinates of quadrature points on initial mesh configuration.
    XQ1, XQ2 = SQ.tabulate_dof_coordinates().reshape((-1, nd)).T  # coordinates
    nQ_local = len(XQ1)

## Stress calculation and degradation function update
Considering the damage in the material, a degradation function need to apply to the stress calculation.

    def sigma(v,deg):
	    strain = eps(v)
	    return (2*mu*strain + lam*tr(strain)*Identity(nd))*deg

Here we define a quadratic degradation function, which is updated before solving the elasticity equation in each iteration based on damage field of the last iteration. 

    ### update degradation function with damage value from previous step
    deg.vector()[:]= np.fmax((1-np_array(dmg))**2,1e-6)

The minimum of the degradation function is a specified small value ($10^{-6}$) instead of zero to prevent numerical instability.

## Strain Energy based History Variable Calculation
In each iteration, the displacement solution can be obtained by solving the elasticity equation with the degraded stress. To calculate the strain energy for the history variable, we first obtain the strain value at each Gaussian integration points in Numpy arrays by

    strain = eps(disp)
    e11= project(strain[0,0],SQ,form_compiler_parameters=form_params).vector().get_local()
    e22= project(strain[1,1],SQ,form_compiler_parameters=form_params).vector().get_local()
    e12= project(strain[0,1],SQ,form_compiler_parameters=form_params).vector().get_local()

One way to calculate the history variable is by iterating over each quad point, solving the principal strain, and computing the corresponding elastic strain energy at all quad points. The iteration way is simple for coding but inefficient, especially for large-scale computing. The vectorized operation can be applied to accelerate the calculation.

We first define two array to store the strain energy value and the strain matrix.

    energy = np.zeros(nQ_local) # energy function
    ### Store the strain tensor on all quad points. (n_quad_points*2*2)
    all_matrix_strain = np.zeros((nQ_local,2,2))
    all_matrix_strain[:,0,0],all_matrix_strain[:,0,1],all_matrix_strain[:,1,0],all_matrix_strain[:,1,1] =\
	    e11,e12,e12,e22

Then the principle strain at each quad points can be obtained by solving the eigenvalues.

    eigs,vecs= np.linalg.eigh(all_matrix_strain)
    eig_1,eig_2 = eigs.max(axis = 1),eigs.min(axis = 1)

 - Miehe's scheme
[Miehe et al.](https://www.sciencedirect.com/science/article/abs/pii/S0045782510001283) defined the positive component of elastic strain energy $\psi_e^+$ by
$$\psi_e^+(\mathbf{\epsilon}) = \frac{\lambda}{2} {\langle \text{tr}(\mathbf{\epsilon}) \rangle}^2 + \mu \text{tr}({\langle \mathbf{\epsilon} \rangle}^2),$$
where $\langle \cdot \rangle = (\ \cdot\ + |\cdot|)/2$ denotes the Macaulay brackets.
- Lo's scheme
[Lo et al.](https://www.sciencedirect.com/science/article/abs/pii/S0022509619306568) proposed a different decomposition scheme in three dimensions. [Sun et al.](https://www.sciencedirect.com/science/article/abs/pii/S2352431621000626) adapted the formulation for the two-dimensional (2-D) plane strain case as
$$
 \begin{cases}
\textbf{if}\ \boldsymbol{\varepsilon}_1 \geq \boldsymbol{\varepsilon}_2 \geq 0,\\
\ \ \textbf{then}\ \psi_e^+ = \dfrac{\lambda}{2}(\boldsymbol{\varepsilon}_1+\boldsymbol{\varepsilon}_2)^2+\mu (\boldsymbol{\varepsilon}_1^2+\boldsymbol{\varepsilon}_2^2);\\
\textbf{elseif}  \ \boldsymbol{\varepsilon}_1 \geq 0 \geq \boldsymbol{\varepsilon}_2 \ \textbf{and} \ (1-\nu) \boldsymbol{\varepsilon}_1 + \nu \boldsymbol{\varepsilon}_2 \geq 0,\\
\ \ \textbf{then}\ \psi_e^+ = \dfrac{E[(1-\nu)\boldsymbol{\varepsilon}_1+\nu \boldsymbol{\varepsilon}_2]^2}{2(1-2\nu)(1-\nu^2)};\\
\textbf{else}\ \psi_e^+ = 0\,
\end{cases} 
$$

Based on the formulation above, we can calculate the strain energy value at each quad point and save them in a Numpy array. Then, we can update the history variable by saving the max value in the loading history.

    ### update history variable
    kappa= np.fmax(energy,kappa)
    phi.vector()[:]= kappa
    
## Solving the Phase Field Damage Evolution Equation

Based on the weak form of the phase field damage evolution equation, the `LHS` and `RHS` can be defined by

    Lhs= ((Constant(eta/delta_t +1/lc)+2*phi)*inner(d,omega)*dx 
          + Constant(lc)*inner(nabla_grad(d),nabla_grad(omega))*dx)
    Rhs= (inner(omega,2*phi)*dx + Constant(eta/delta_t)*inner(omega,dmg)*dx)

Then, the equation can be solved by

    ### solve the variational problem
    w = Function(S)
    solve(Lhs==Rhs,w,form_compiler_parameters=form_params)
The damage field defined on element nodes is stored in `w`. To get the damage field defined on the quad points, a interpolate method of FEniCS is applied. We also need to ensure that the damage value after interpolation should in $[0,1]$ and irreversible.

    dmgv = np.clip(interpolate(w, SQ).vector().get_local(), 0, 1) # Interpolate damage results from nodes to quad points
    dmg.vector()[:] = np.fmax(dmgv, dmg.vector().get_local()) # Damage evolution is irreversible

Then the damage field can be used for the next iteration calculation.
## Loading speed
The displacement increment should be sufficiently small for obtaining accurate results, especially after the crack starts to propagate. Please refer to Appendix B in [Sun et al.](https://www.sciencedirect.com/science/article/abs/pii/S2352431621000626)
## Stopping Criteria

We can track the reaction forces at each iteration. When it reach to the maximum and drop to a very small value, the loop can be stopped to prevent unnecessary calculation.

## Results
![The damage field](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/Damage_PhaseField.jpg?raw=true)
![The load-displacement curve](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/Load_Disp_PhaseField.jpg?raw=true)

