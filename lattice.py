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
  
def update(grid):
  grid=move(grid)
  #grid=pressure(grid) ???
  #grid=relax(grid)
  return grid
  