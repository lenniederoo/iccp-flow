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
plt.figure(1)
plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]),'r+')
plt.xlabel('Velocity in x direction') 
plt.ylabel('Position in y') 
plt.title('Flow profile in pipe')
plt.show()