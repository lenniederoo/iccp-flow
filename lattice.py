import numpy as np
import copy


def init_grid(Nxgrid,Nygrid,dens):
  grid=np.ones((Nxgrid,Nygrid,9),dtype=float)*dens/9
  return grid
  
def move(grid):
  print 'moving'
  clonedgrid=copy.deepcopy(grid)
  grid=np.zeros(clonedgrid.shape,dtype=float)
  for j in xrange(0,grid.shape[0]):
    for k in xrange(0,grid.shape[1]):    
        if (0<k<grid.shape[1]-1):
          grid[j,k,0]+=clonedgrid[j,k,0]
          grid[j,k,1]+=clonedgrid[(j-1)%(grid.shape[0]-1),k,1]
          grid[j,k,2]+=clonedgrid[(j-1)%(grid.shape[0]-1),k-1,2]
          grid[j,k,3]+=clonedgrid[j,k-1,3]
          grid[j,k,4]+=clonedgrid[(j+1)%(grid.shape[0]-1),k-1,4]
          grid[j,k,5]+=clonedgrid[(j+1)%(grid.shape[0]-1),k,5]
          grid[j,k,6]+=clonedgrid[(j+1)%(grid.shape[0]-1),k+1,6]
          grid[j,k,7]+=clonedgrid[j,k+1,7]
          grid[j,k,8]+=clonedgrid[(j-1)%(grid.shape[0]-1),k+1,8]
        elif k==0:
          grid[j,k,0]+=clonedgrid[j,k,0]
          grid[j,k,1]+=clonedgrid[(j-1)%(grid.shape[0]-1),k,1]
          grid[j,k,2]+=clonedgrid[(j-1)%(grid.shape[0]-1),k,6]
          grid[j,k,3]+=clonedgrid[j,k,7]
          grid[j,k,4]+=clonedgrid[(j+1)%(grid.shape[0]-1),k,8]
          grid[j,k,5]+=clonedgrid[(j+1)%(grid.shape[0]-1),k,5]
          grid[j,k,6]+=clonedgrid[(j+1)%(grid.shape[0]-1),k+1,6]
          grid[j,k,7]+=clonedgrid[j,k+1,7]
          grid[j,k,8]+=clonedgrid[(j-1)%(grid.shape[0]-1),k+1,8]
        else:
          grid[j,k,0]+=clonedgrid[j,k,0]
          grid[j,k,1]+=clonedgrid[(j-1)%(grid.shape[0]-1),k,1]
          grid[j,k,2]+=clonedgrid[(j-1)%(grid.shape[0]-1),k-1,2]
          grid[j,k,3]+=clonedgrid[j,k-1,3]
          grid[j,k,4]+=clonedgrid[(j+1)%(grid.shape[0]-1),k-1,4]
          grid[j,k,5]+=clonedgrid[(j+1)%(grid.shape[0]-1),k,5]
          grid[j,k,6]+=clonedgrid[(j+1)%(grid.shape[0]-1),k,2]
          grid[j,k,7]+=clonedgrid[j,k,3]
          grid[j,k,8]+=clonedgrid[(j-1)%(grid.shape[0]-1),k,4]
  return grid
  
def update(grid,relaxt,pressgradvel):
  grid=move(grid)  
  grid=relax(grid,relaxt,pressgradvel)
  return grid
  
def calc_velocity(grid):
    velocity=np.zeros((grid.shape[0],grid.shape[1],2),dtype=float)
    velocity[:,:,0]=grid[:,:,1]+0.5*np.sqrt(2)*(grid[:,:,2]+grid[:,:,8]-grid[:,:,4]-grid[:,:,6])-grid[:,:,5]    #x-direction
    velocity[:,:,1]=grid[:,:,3]-grid[:,:,7]+0.5*np.sqrt(2)*(grid[:,:,2]+grid[:,:,4]-grid[:,:,6]-grid[:,:,8])    #y-direction
    return velocity
    
def add_pressure_grad(grid,pressgradvel):
    velocity=calc_velocity(grid)
    velocity[:,:,0]+=pressgradvel
    return velocity
    
def calc_eq(grid,pressgradvel):
  velocity=add_pressure_grad(grid,pressgradvel)
  grideq=grid  
  return grideq
  
def relax(grid,relaxt,pressgradvel):
    grideq=calc_eq(grid,pressgradvel)
    grid=(1-(1/relaxt))*grid+grideq*(1/relaxt)
    return grid
    