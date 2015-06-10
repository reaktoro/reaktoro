import numpy as npy
from dolfin import *


class Transport:
    def __init__(self, u, **kwargs):
        self.dt = Constant(0.0)

        velocity  = kwargs.get("velocity", Constant(0.0))
        diffusion = kwargs.get("diffusion", Constant(0.0))
        source = kwargs.get("source", Constant(0.0))

        V = u.function_space()
        mesh = V.mesh()

        u0 = u

        u = TrialFunction(V)
        v = TestFunction(V)

        dt = self.dt

        # Mid-point solution
        u_mid = 0.5*(u0 + u)

        # Residual
        r = u - u0 + dt*(dot(velocity, grad(u_mid)) - \
            diffusion*div(grad(u_mid)) - source)

        # Galerkin variational problem
        F = v*(u - u0)*dx + dt*(v*div(u_mid*velocity)*dx + \
            diffusion*dot(grad(v), grad(u_mid))*dx - source*v*dx)

        # Add SUPG stabilisation terms
        h = CellSize(mesh)
        vnorm = sqrt(dot(velocity, velocity))
        tau = pow((2.0/dt)**2 + (2.0*vnorm/h)**2 + 9*(4*diffusion/h**2), -0.5)

        F += tau*dot(velocity, grad(v))*r*dx

        # Create bilinear and linear forms
        self.a = lhs(F)
        self.L = rhs(F)

        self.A = dolfin.PETScMatrix()
        self.b = dolfin.PETScVector()

        self.solver = dolfin.LUSolver()
        self.solver.parameters["same_nonzero_pattern"] = True

    def step(self, u, dt, bc):
        self.dt.assign(dt)
        assemble(self.a, tensor=self.A)
        assemble(self.L, tensor=self.b)
        bc.apply(self.A)
        bc.apply(self.b)
        self.solver.set_operator(self.A)
        self.solver.solve(u.vector(), self.b)



#
# class Transport:
#     def __init__(self, mesh):
#         self.mesh = mesh
#         self.V = FunctionSpace(mesh, "CG", 1)
#         self.Q = VectorFunctionSpace(mesh, "CG", 2)
#
#         self.velocity = Function(self.Q)
#         self.diffusion = Function(self.V)
#         self.source = Function(self.V)
#
#         u = TrialFunction(self.V)
#         v = TestFunction(self.V)
#
#         self.u0 = Function(self.V)
#         self.dt = Constant(0.0)
#
#         # Mid-point solution
#         u_mid = 0.5*(self.u0 + u)
#
#         # Residual
#         r = u - self.u0 + \
#             self.dt*(dot(self.velocity, grad(u_mid)) -
#                 self.diffusion*div(grad(u_mid)) - self.source)
#
#         # Galerkin variational problem
#         F = v*(u - self.u0)*dx + \
#             self.dt*(v*div(u_mid*self.velocity)*dx +
#                 self.diffusion*dot(grad(v), grad(u_mid))*dx -
#                     self.source*v*dx)
#
#         # Add SUPG stabilisation terms
#         h = CellSize(mesh)
#         vnorm = sqrt(dot(self.velocity, self.velocity))
#         tau = pow((2.0/self.dt)**2 + (2.0*vnorm/h)**2 + 9*(4*self.diffusion/h**2), -0.5)
#
#         F += tau*dot(self.velocity, grad(v))*r*dx
#
#         # Create bilinear and linear forms
#         self.a = lhs(F)
#         self.L = rhs(F)
#
#         self.A = dolfin.PETScMatrix()
#         self.b = dolfin.PETScVector()
#
#         self.solver = dolfin.LUSolver()
#         self.solver.parameters["same_nonzero_pattern"] = True
#
#     def setInitialCondition(self, u0):
#         self.u0.assign(u0)
#
#     def setInitialConditionAsArray(self, u0):
#         self.u0.vector()[:] = npy.array(u0)
#
#     def setVelocity(self, velocity):
#         if type(velocity) not in [Constant, Function]:
#             velocity = project(velocity, self.Q)
#         self.velocity.assign(velocity)
#         self.updated = True
#
#     def setDiffusionCoeff(self, diffusion):
#         if type(diffusion) not in [Constant, Function]:
#             diffusion = project(diffusion, self.V)
#         self.diffusion.assign(diffusion)
#         self.updated = True
#
#     def setSource(self, source):
#         self.source.assign(source)
#         self.updated = True
#
#     def setBoundaryCondition(self, bc):
#         self.bc = bc
#
#     def u(self):
#         return self.u0
#
#     def step(self, dt):
#         self.dt.assign(dt)
#
#         if self.updated:
#             assemble(self.a, tensor=self.A)
#             self.bc.apply(self.A)
#             self.solver.set_operator(self.A)
#             self.updated = False
#
#         assemble(self.L, tensor=self.b)
#
#         self.bc.apply(self.b)
#
#         self.solver.solve(self.u0.vector(), self.b)
#         self.solver.parameters['reuse_factorization'] = True
