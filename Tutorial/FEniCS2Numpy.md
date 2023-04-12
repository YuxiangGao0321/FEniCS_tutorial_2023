# Save FEniCS results as a Numpy array
We consider a displacement field $u$ defined on a vector function space `V` where

    V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)
    u = Function(V)
## Results at nodes (Displacement result in each component)
To extract the displacement field at each component, we need to define a scaler function space `S` on the nodes by

    S = FunctionSpace(mesh,"Lagrange",1)

and project each component of the vector field to the scaler field by

    ux = project(u[0],S)

Then we can convert the FEniCS function variable to a Numpy array by

    ux_array = ux.vector().get_local()
The x-direction displacement results at each node are saved in the Numpy array `ux_array`. The length of `ux_array` is the number of the node.
## Results at integration points (Strain/Stress result in each component)
As we mentioned before, in FEniCS, the strain and stress tensor can be calculated by

    strain = sym(grad(u))
    stress = lam*tr(strain)*Identity(2) + 2.0*mu*strain

To extract the strain/stress result in each component, we need to define a scaler function space `SQ` on the integration points by

    deg_quad = 1
	scalar_quad = FiniteElement("Quadrature", cell=mesh.ufl_cell(),
                            degree=deg_quad, quad_scheme="default")
	SQ = FunctionSpace(mesh, scalar_quad)  # quadrature points in scalar space
	form_params = {"quadrature_degree": deg_quad}

and project each component of the strain/stress tensor to the scaler field on the integration points by

    sxx = project(stress[0,0],SQ,form_compiler_parameters=form_params)

Then we can convert the FEniCS function variable to a Numpy array by

    sxx_array = sxx.vector().get_local()
The $\sigma_{xx}$ at each integration points are saved in the Numpy array `x_array`. The length of `x_array` is the number of the integration points.

## Save and load the results

The Numpy array can be saved as a text file by [numpy.savetxt](https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html). The text field can be load by [numpy.loadtxt](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html)

