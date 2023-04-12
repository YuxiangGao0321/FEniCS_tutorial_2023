from dolfin import *
from matplotlib import pyplot as plt
from fenics import *
# Load mesh and subdomains
# mesh = Mesh("dolfin_fine.xml.gz")
# sub_domains = MeshFunction("size_t", mesh, "dolfin_fine_subdomains.xml.gz")
# plot(mesh)
# plot(sub_domains)
i = 10
mesh = UnitSquareMesh(i, i)
# Define function spaces
# use the MixedElement class to define a mixed element and then pass this to the FunctionSpace constructor to create a mixed function space
# Define function spaces
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)
#"CG" for continuous Galerkin elements), the cell type (e.g., mesh.ufl_cell() for the cell type of the mesh- triangles ), and the polynomial degree of the finite element space.

# V = VectorFunctionSpace(mesh, "CG", 2)
# Q = FunctionSpace(mesh, "CG", 1)
# W = V * Q
# noslip = Constant((0, 0))
noslip = project(Constant((0, 0)), W.sub(0).collapse())
bcs = DirichletBC(W.sub(0), noslip, "on_boundary")
# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
# f = Constant((0, 0))
# Define source term
f = project(Expression(("-(4*x[1]*(1-x[1])*(2*x[1]-1)*((1-2*x[0])*(1-2*x[0])-2*x[0]*(1-x[0])) + \
                  12*x[0]*x[0]*(1-x[0])*(1-x[0])*(1-2*x[1])) + (1-2*x[0])*(1-x[1])",
                "-(4*x[0]*(1-x[0])*(1-2*x[0])*((1-2*x[1])*(1-2*x[1])-2*x[1]*(1-x[1])) + \
                  12*x[1]*x[1]*(1-x[1])*(1-x[1])*(2*x[0]-1)) - x[0]*(1-x[0])"), degree=2),W.sub(0).collapse())
a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = inner(f, v)*dx
# Compute solution
w = Function(W)
solve(a == L, w, bcs)
# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)
print "Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2")
print "Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2")
# # Split the mixed solution using a shallow copy
(u, p) = w.split()
# Save solution in VTK format
ufile_pvd = File("velocityTH.pvd")
ufile_pvd << u
pfile_pvd = File("pressureTH.pvd")
pfile_pvd << p
# Compute L2 error with respect to analytical solutions
ua = project(Expression(("x[0]*x[0]*(1-x[0])*(1-x[0])*2*x[1]*(1-x[1])*(2*x[1]-1)",
                 "x[1]*x[1]*(1-x[1])*(1-x[1])*2*x[0]*(1-x[0])*(1-2*x[0])"),
                degree=3),W.sub(0).collapse())

pa = project(Expression("x[0]*(1-x[0])*(1-x[1])-0.0833333", degree=1),W.sub(1).collapse())
mean_p = assemble(p*dx(domain=mesh)) / assemble(1*dx(domain=mesh))
print("Mean pressure:", mean_p)
#  mean pressure from p
mean_p = project(mean_p,W.sub(1).collapse())
p_minus_mean = p - mean_p
p_minus_mean = project(p_minus_mean,W.sub(1).collapse())
# compute erors
err_u_L2 = errornorm(u, ua, "L2")
err_p_L2 = errornorm(p_minus_mean, pa, "L2")

err_u_H1 = errornorm(u, ua, "H10")
err_p_H1 = errornorm(p_minus_mean, pa, "H10")

# Compute L2 error with respect to analytical solutions

print("L2 error in velocity: %.9f" % err_u_L2)
print("L2 error in pressure: %.9f" % err_p_L2)
print("H1 error in velocity: %.9f" % err_u_H1)
print("H1 error in pressure: %.9f" % err_p_H1)



# Plot solution
u_plot = plot(u)
plt.colorbar(u_plot)
plt.title("velocity - Mini elements")
plt.savefig("uTH.jpg",dpi = 500)
plt.close()

p_plot = plot(p_minus_mean)
plt.colorbar(p_plot)
plt.title("p")
plt.savefig("pTH.jpg",dpi = 500)
plt.close()


ua_plot = plot(ua)
plt.colorbar(ua_plot)
plt.title("ua")
plt.savefig("uaTH.jpg",dpi = 500)
plt.close()


pa_plot = plot(pa)
plt.colorbar(pa_plot)
plt.title("pa")
plt.savefig("paTH.jpg",dpi = 500)
plt.close()