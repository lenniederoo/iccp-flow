from lattice import *
print 'hello world'

Nxgrid=10 #gridpoints x direction
Nygrid=4  #gridpoints y direction
dens=900
timesteps=10
pressgradvel=0.2
relaxt=0.5

grid=init_grid(Nxgrid,Nygrid,dens)
for i in xrange(0,timesteps):
  grid=update(grid,relaxt,pressgradvel)
  
print grid