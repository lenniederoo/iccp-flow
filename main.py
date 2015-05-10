from lattice import *
import numpy as np
from pylab import *

Nxgrid=100 #gridpoints x direction
Nygrid=50  #gridpoints y direction
dens=1
timesteps=1000
pressgradvel=0.005
relaxt=3.3
blocks = np.array([[30,30],[60,40],[80,20]]) # Adds 1 by 1 blocks on position [x,y]
bigblocks = np.array([[10,20],[20,30]]) # Adds bigger blocks of size [xstart,xend] by [ystart,yend]
e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]])

grid=init_grid(Nxgrid,Nygrid,dens,blocks,bigblocks)
for i in xrange(0,timesteps):
    grid=update(grid,relaxt,pressgradvel,blocks,bigblocks,e)
#print grid
velocity=calc_velocity(grid)
plt.figure(1)
plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]),'r+')
plt.xlabel('Velocity in x direction') 
plt.ylabel('Position in y') 
plt.title('Flow profile in pipe')
plt.show()

plt.figure(2)
tempvel = np.copy(velocity)
v = np.transpose(tempvel)
Q = quiver(v[0,1:-1,:], v[1,1:-1,:])
l,r,b,t = axis()
dx, dy = r-l, t-b
axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])
title('velocity profile')