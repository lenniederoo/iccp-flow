from lattice import *
import numpy as np
#print 'hello world'

Nxgrid=100 #gridpoints x direction
Nygrid=50  #gridpoints y direction
dens=1
timesteps=1000
pressgradvel=0.005
relaxt=1.85
blocks = np.array([[40,20]]) # Adds 1 by 1 blocks on position [x,y]


grid=init_grid(Nxgrid,Nygrid,dens,blocks)
for i in xrange(0,timesteps):
    grid=update(grid,relaxt,pressgradvel,blocks)
#print grid
velocity=calc_velocity(grid)
plt.figure(1)
plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]),'r+')
plt.xlabel('Velocity in x direction') 
plt.ylabel('Position in y') 
plt.title('Flow profile in pipe')
plt.show()