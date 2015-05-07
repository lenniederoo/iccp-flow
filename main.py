from lattice import *
#print 'hello world'

Nxgrid=100 #gridpoints x direction
Nygrid=50  #gridpoints y direction
dens=1
timesteps=1000
pressgradvel=0.005
relaxt=1.85

grid=init_grid(Nxgrid,Nygrid,dens)
for i in xrange(0,timesteps):
    grid=update(grid,relaxt,pressgradvel)
#print grid
velocity=calc_velocity(grid)
plt.figure(2)
plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]),'r+')
plt.show()