# Lattice movement
import numpy as np
import matplotlib.pyplot as plt
import copy

def init_grid(Nxgrid,Nygrid,dens,blocks):
  grid=np.ones((Nxgrid,Nygrid,9),dtype=float)*dens/9
  for i in range (len(blocks)):
    grid[blocks[i,0],blocks[i,1],:]=0
  return grid

def move(grid,e):
  for i in range(9):
    grid[:,:,i]=np.roll(np.roll(grid[:,:,i],e[0,i],axis=0),e[1,i],axis=1)
  clonedgrid=copy.deepcopy(grid)
  for i in range(9):
    grid[:,0,i]=clonedgrid[:,0,((i+3)%8) +1]
    grid[:,grid.shape[1]-1,i]=clonedgrid[:,grid.shape[1]-1,((i+3)%8) +1]
  return grid  
  
def update(grid,relaxt,pressgradvel,blocks,e):
  grid=move(grid,e)  
  grid=relax(grid,relaxt,pressgradvel,blocks,e)
  return grid
  
def calc_velocity(grid):
    velocity=np.zeros((grid.shape[0],grid.shape[1],2),dtype=float)
    velocity[:,:,0]=(1/np.sum(grid,axis=2))*(grid[:,:,1]+(grid[:,:,2]+grid[:,:,8]-grid[:,:,4]-grid[:,:,6])-grid[:,:,5])    #x-direction
    velocity[:,:,1]=(1/np.sum(grid,axis=2))*(grid[:,:,3]-grid[:,:,7]+(grid[:,:,2]+grid[:,:,4]-grid[:,:,6]-grid[:,:,8]))    #y-direction
    return velocity
    
def add_pressure_grad(grid,pressgradvel,blocks):
    velocity=calc_velocity(grid)
    velocity[:,1:velocity.shape[1]-1,0]+=pressgradvel
    for p in range (len(blocks)):  
      velocity[blocks[p,0],blocks[p,1],0]-=pressgradvel
    return velocity
    
def calc_eq(grid,pressgradvel,blocks,e):
  velocity=add_pressure_grad(grid,pressgradvel,blocks)
  Velsqr=velocity[:,:,0]**2+velocity[:,:,1]**2
  grideq=np.zeros(grid.shape,dtype=float)
  w= np.array([4./9, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36])
  Prod=np.dot(e.T,velocity.transpose(0,2,1))
  for i in xrange(0,9):
    grideq[:,:,i]=w[i]*np.sum(grid,axis=2)*(1+3* Prod[i] + 9/2.0* Prod[i]**2 -3.0/2 * Velsqr)
  return grideq
  
def relax(grid,relaxt,pressgradvel,blocks,e):
    grideq=calc_eq(grid,pressgradvel,blocks,e)
    grid=(1-(1/relaxt))*grid+grideq*(1/relaxt)
    return grid
    