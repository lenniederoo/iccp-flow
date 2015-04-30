# Lattice movement
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
          grid[j,k,0]+=0
          grid[j,k,1]+=0
          grid[j,k,2]+=clonedgrid[(j+1)%(grid.shape[0]-1),k+1,6]
          grid[j,k,3]+=clonedgrid[j,k+1,7]
          grid[j,k,4]+=clonedgrid[(j-1)%(grid.shape[0]-1),k+1,8]
          grid[j,k,5]+=0
          grid[j,k,6]+=0
          grid[j,k,7]+=0
          grid[j,k,8]+=0
        else:
          grid[j,k,0]+=0
          grid[j,k,1]+=0
          grid[j,k,2]+=0
          grid[j,k,3]+=0
          grid[j,k,4]+=0
          grid[j,k,5]+=0
          grid[j,k,6]+=clonedgrid[(j-1)%(grid.shape[0]-1),k-1,2]
          grid[j,k,7]+=clonedgrid[j,k-1,3]
          grid[j,k,8]+=clonedgrid[(j+1)%(grid.shape[0]-1),k-1,4]
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
  grideq=np.zeros(grid.shape,dtype=float)
  for i in xrange(0,9):
      if i==0:
          w=4/9.0
          ea=0
          eb=0
      elif i==2:
           w=1.0/9
           ea=1
           eb=1
      elif i==4:
          w=1.0/9
          ea=-1
          eb=1
      elif i==6:
          w=1.0/9
          ea=-1
          eb=-1
      elif i==8:
          w=1.0/9
          ea=1
          eb=-1
      elif i==1:
          w=1.0/36
          ea=1
          eb=0
      elif i==3:
          w=1.0/36
          ea=0
          eb=1
      elif i==5:
          w=1.0/36
          ea=-1
          eb=0
      elif i==7:
          w=1.0/36
          ea=0
          eb=-1
      grideq[:,:,i]=w*np.sum(grid,axis=2)*(1+3*ea*velocity[:,:,0]+9/2.0*ea*eb*velocity[:,:,0]*velocity[:,:,1]-3.0/2*velocity[:,:,0]*velocity[:,:,0])
  return grideq
  
def relax(grid,relaxt,pressgradvel):
    grideq=calc_eq(grid,pressgradvel)
    grid=(1-(1/relaxt))*grid+grideq*(1/relaxt)
    return grid
    