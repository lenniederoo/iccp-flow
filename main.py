from lattice import *
print 'hello world'

Nxgrid=100 #gridpoints x direction
Nygrid=10  #gridpoints y direction
dens=900
timesteps=100
pressgradvel=0.2
relaxt=0.5

grid=init_grid(Nxgrid,Nygrid,dens)
for i in xrange(0,timesteps):
  grid=update(grid,relaxt,pressgradvel)
print grid