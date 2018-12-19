import reaktoro as rkt
import matplotlib.pyplot as plt
import numpy as np

# Auxiliary time related constant
day = 86400 #[s]

# Parameters for the transport simulation
nsteps = 225      # the number of steps in the transport simulation
ncells = 100      # the number of cells in the discretization
xl = 0.0          # the x-coordinate of the left boundary
xr = 100.0        # the x-coordinate of the right boundary
D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the velocity (in units of m/s)
dt = 0.5*day      # the time step (in units of s)
ul = 1            # concentration at the left boundary (mol/m3)

mesh = rkt.Mesh(ncells, xl, xr)

x = mesh.xcells()
 
transport = rkt.TransportSolver()
transport.setMesh(mesh)
transport.setVelocity(v)
transport.setDiffusionCoeff(D)
transport.setBoundaryValue(ul)
transport.setTimeStep(dt)
      
transport.initialize()
  
u = np.zeros(ncells)

useries = []
tseries = []

for i in range(1,nsteps):

    # Record at every 25 steps the current time and concentration
    if i % 25 == 0:
        tseries.append(i * dt / day)
        useries.append(u.copy())

    # Perform one time step
    transport.step(u)

fig, ax = plt.subplots(figsize=(10, 20))    
ax.set(xlabel='x [m]', ylabel='u [mol/m3]')

for i in range(len(tseries)):
    ax.plot(x, useries[i], label='{} days'.format(tseries[i]))

legend = ax.legend(fontsize='x-large')

fig.set_size_inches((18, 7))

plt.show()
