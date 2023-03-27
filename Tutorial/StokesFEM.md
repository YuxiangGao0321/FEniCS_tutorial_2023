# Some methods for solving Stokes equation with FEM
Various mixed finite element pairs have been proposed for solving the Stokes problem and other problems in numerical fluid dynamics. These pairs differ in their global regularity and the local order of the polynomials they use. Mixed finite element pairs are typically denoted by the symbol $[\mathcal{P}_k]^d/\mathcal{P}_m$, which indicates the number of space dimensions $d$ (in this case, 2) and the local order of the polynomials used for the velocity ($k$) and pressure ($m$) spaces. As is well known, the pair $[\mathcal{P}_1]^2/\mathcal{P}_1$ that would be computationally efficient is however unstable.
## MINI element
The MINI finite element was introduced by Arnold, Brezzi, and Fortin as a specific discretization method for the Stokes problem. It is characterized by the $[\mathcal{P}_{1+b}]^2/\mathcal{P}_1$ mixed finite element pair, in which the discrete velocity space $X_h$ is comprised of continuous, piecewise-linear polynomials enriched with local cubic bubble functions, while the discrete pressure space $M_h$ is comprised of continuous, piecewise-linear polynomials.

The MINI element was developed as a means of stabilizing the $[\mathcal{P}_1]^2/\mathcal{P}_1$ mixed finite element pair, which is known to be unstable. The addition of the cubic bubble functions to the discrete velocity space serves to enrich it, thus achieving greater stability while maintaining computational efficiency. The bubble function used in the MINI element is a cubic polynomial that is defined locally in each triangle and is given by the product of the barycentric coordinates $\varphi_1, \varphi_2,$ and $\varphi_3$ of the triangle itself. The form of the MINI element is shown below.

$$
\begin{align}
v_{h} = a_h^{(v)}  + b_h^{(v)}x + c_h^{(v)}y + d_h^{(v)} \varphi_1(x,y) \varphi_2(x,y) \varphi_3(x,y) \\
q_{h} = a_h^{(q)} + b_h^{(q)} x + c_h^{(q)} y
\end{align}$$

Where $a_h, b_h, c_h, d_h$ are unknown constants.
![MINI elements](https://github.com/YuxiangGao0321/FEniCS_tutorial_2023/blob/main/figs/MNIelement.jpg?raw=true|width=80)

To define the `FunctionSpace` of the MINI element `(P1+B)*Q` in FEniCS 2016.2.0, one can follow the [documentation](https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/demo/documented/stokes-mini/python/documentation.html) and use the following code:

    ### Build function spaces on Mini element
    P1 = VectorElement("Lagrange", mesh.ufl_cell(), 1)
    B = VectorElement("Bubble",   mesh.ufl_cell(), 3)
    Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, (P1 + B) * Q)

The bilinear and linear forms corresponding to the weak mixed formulation of the Stokes equations are defined as follows:

    ### Define variational problem
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    f = Constant((0, 0))
    a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
    L = inner(f, v)*dx

## Stabilized first order elements

As previously noted, using the computationally efficient pair $[\mathcal{P}_1]^2/\mathcal{P}_1$ is not viable due to instability issues. These first-order elements are susceptible to spurious pressure oscillations in regions of high vorticity or sharp velocity gradients.

To address this issue, stabilization techniques can be employed in the finite element method. In stabilized first-order elements, a stabilizing term is added to the momentum equation, which incorporates the gradient of the velocity and pressure fields, as well as a stabilization parameter. The value of the parameter can be adjusted to achieve the desired level of stabilization.

In FEniCS 2016.2.0, we can first define the $[\mathcal{P}_1]^2/\mathcal{P}_1$ `FunctionSpace` by 

    ### Define function spaces
    P1_v = VectorElement("Lagrange", mesh.ufl_cell(), 1)
    P1_s = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    P1P1 = P1_v * P1_s
    system = FunctionSpace(mesh, P1P1)

The bilinear and linear forms corresponding to the stabilized weak mixed formulation of the Stokes equations are defined as follows:
    
    ### Define variational problem
    (v, q) = TestFunctions(system)
    (u, p) = TrialFunctions(system)
    f = Constant((0, 0))
    h = CellSize(mesh)
    beta  = 0.2
    delta = beta*h*h
    a = (inner(grad(v), grad(u)) - div(v)*p + q*div(u) + \
        delta*inner(grad(q), grad(p)))*dx
    L = inner(v + delta*grad(q), f)*dx

## Taylor-Hood elements
The Taylor-Hood elements (also known as $[\mathcal{P}_{2}]^2/\mathcal{P}_1$ element) use piecewise quadratic polynomials for the velocity field and piecewise linear polynomials for the pressure field. This means that the velocity is approximated by a quadratic polynomial within each element of the mesh, while the pressure is approximated by a linear polynomial within each element.

To define a `FunctionSpace` for the Taylor-Hood elements in FEniCS 2016.2.0, one can use the following code, as specified in the [documentation](https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/demo/documented/stokes-taylor-hood/python/documentation.html):

    ### Define function spaces
    P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    TH = P2 * P1
    W = FunctionSpace(mesh, TH)

The weak form is defined in the same way as the MINI element.
### Reference

- Cioncolini, A., & Boffi, D. (2019). The MINI mixed finite element for the Stokes problem: An experimental investigation. _Computers & Mathematics with Applications_, _77_(9), 2432–2446. [https://doi.org/10.1016/j.camwa.2018.12.028](https://doi.org/10.1016/j.camwa.2018.12.028)
- Anjos, G. R. (2020). Numerical Investigation of Two-Phase Flows in Corrugated Channel with Single and Multiples Drops. _Fluids_, _6_(1), 13. [https://doi.org/10.3390/fluids6010013](https://doi.org/10.3390/fluids6010013)


