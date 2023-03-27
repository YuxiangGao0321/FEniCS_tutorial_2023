from dolfin import *
from matplotlib import pyplot as plt
from fenics import *

# Geom
L = 25.
H = 1.
Nx = 250
Ny = 10
mesh = RectangleMesh(Point(0., 0.), Point(L, H), Nx, Ny, "crossed")
# mesh = Mesh("mesh.xml")

# plot(mesh, lw = 0.1)
# plt.savefig('mesh.jpg',dpi = 500)

# material parameter
E = Constant(1e5)
nu = Constant(0.3)
model = "plane_stress"

mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)
if model == "plane_stress":
    lmbda = 2*mu*lmbda/(lmbda+2*mu)

# Strain
def eps(v):
    return sym(grad(v))

# Stress
def sigma(v):
    return lmbda*tr(eps(v))*Identity(2) + 2.0*mu*eps(v)

# Body force vector
rho_g = 1e-3
f = Constant((0, -rho_g))


# Neumann BC
# class right_edge(SubDomain):
#     """Boundary on the right domain edge."""
#     def inside(self, x, on_boundary):
#         return near(x[0], L) and on_boundary

# right= right_edge()
# boundaries = FacetFunction("size_t", mesh, 0)

# right.mark(boundaries, 1)

# ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# Define function space
V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
du = TrialFunction(V)
u_ = TestFunction(V)

# LHS and RHS of weak form
a = inner(sigma(du), eps(u_))*dx
l = inner(f, u_)*dx

# a = inner(sigma(du), eps(u_))*dx
# Tx = 10
# l = inner(f, u_)*dx + inner(u_,Constant((Tx,0)))*ds(1)

# fix left boundary
def left(x, on_boundary):
    return near(x[0], 0.)

bc = DirichletBC(V, Constant((0.,0.)), left)

# Define a variable to store displacement results
u = Function(V, name="Displacement")
# solve the result
solve(a == l, u, bc)

plot(1e3*u, mode="displacement")
plt.savefig('displacement.jpg',dpi = 300)
plt.close()

# Strain/Stress results
# """Quadrature elements and function spaces."""
# deg_quad = 1
# scalar_quad = FiniteElement("Quadrature", cell=mesh.ufl_cell(),
#                             degree=deg_quad, quad_scheme="default")
# SQ = FunctionSpace(mesh, scalar_quad)  # quadrature points in scalar space
# form_params = {"quadrature_degree": deg_quad}
# #extract strain and stress tensor
# strain = eps(u)
# stress = sigma(u)

# Plot strain
# e11_IP= project(strain[0,0],SQ,form_compiler_parameters=form_params)
# plot(e11_IP)
# plt.savefig('e11.jpg',dpi = 300)

# Convert the results to a numpy array
# e11_IP= project(strain[0,0],SQ,form_compiler_parameters=form_params).vector().get_local()
# e22_IP= project(strain[1,1],SQ,form_compiler_parameters=form_params).vector().get_local()
# e12_IP= project(strain[0,1],SQ,form_compiler_parameters=form_params).vector().get_local()

# s11_IP=project(stress[0,0],SQ,form_compiler_parameters=form_params).vector().get_local()
# s22_IP=project(stress[1,1],SQ,form_compiler_parameters=form_params).vector().get_local()
# s12_IP=project(stress[0,1],SQ,form_compiler_parameters=form_params).vector().get_local()