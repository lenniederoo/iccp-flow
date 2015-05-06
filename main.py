from lattice import *
print 'hello world'

Nxgrid=100 #gridpoints x direction
Nygrid=50  #gridpoints y direction
dens=900
timesteps=100
pressgradvel=0.2
relaxt=1.85

grid=init_grid(Nxgrid,Nygrid,dens)
for i in xrange(0,timesteps):
    grid=update(grid,relaxt,pressgradvel)
print grid
velocity=calc_velocity(grid)
plt.figure(2)
plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]))
plt.show()