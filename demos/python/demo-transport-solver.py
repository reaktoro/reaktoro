import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import reaktoro as rk

# Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Parameters for the transport simulation
nsteps = 225      # the number of steps in the transport simulation
ncells = 100      # the number of cells in the discretization
xl = 0.0          # the x-coordinate of the left boundary
xr = 100.0        # the x-coordinate of the right boundary
D  = 1.0e-9       # the diffusion coefficinet (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 0.5*day      # the time step (in units of s)
ul = 1            # amount at the left boundary

mesh = rk.Mesh(ncells, xl, xr)

x = mesh.xcells()
 
ts = rk.TransportSolver()
ts.setMesh(mesh)
ts.setVelocity(v)
ts.setDiffusionCoeff(D)
ts.setBoundaryValue(ul)
ts.setTimeStep(dt)
      
ts.initialize()
  
numerical_u = np.zeros(ncells)

amounts = []
times = []

for i in range(1,nsteps):
    if (i%25 == 0):
        amounts.append(copy.deepcopy(numerical_u))
        times.append(i*dt/day)
    ts.step(numerical_u)

fig, ax = plt.subplots(figsize=(10, 20))    
ax.set(xlabel='x [m]', ylabel='amount [mol]')

for i in range(0,len(times)):
    ax.plot(x, amounts[i], label='{} days'.format(times[i]))

legend = ax.legend(fontsize='x-large')

fig.set_size_inches((18, 7))

plt.show()
