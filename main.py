from lattice import *
print 'hello world'

Nxgrid=10 #gridpoints x direction
Nygrid=4  #gridpoints y direction
dens=900
timesteps=10

grid=init_grid(Nxgrid,Nygrid,dens)
grid=move(grid)
for i in xrange(0,timesteps):
  update(grid)
print grid