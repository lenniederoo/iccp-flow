# Lattice movement
import numpy as np
import copy
import matplotlib.pyplot as plt


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
          grid[j,k,1]+=clonedgrid[(j-1)%(grid.shape[0]),k,1]
          grid[j,k,2]+=clonedgrid[(j-1)%(grid.shape[0]),k-1,2]
          grid[j,k,3]+=clonedgrid[j,k-1,3]
          grid[j,k,4]+=clonedgrid[(j+1)%(grid.shape[0]),k-1,4]
          grid[j,k,5]+=clonedgrid[(j+1)%(grid.shape[0]),k,5]
          grid[j,k,6]+=clonedgrid[(j+1)%(grid.shape[0]),k+1,6]
          grid[j,k,7]+=clonedgrid[j,k+1,7]
          grid[j,k,8]+=clonedgrid[(j-1)%(grid.shape[0]),k+1,8]
        elif k==0:
          grid[j,k,0]+=0
          grid[j,k,1]+=0
          grid[j,k,2]+=clonedgrid[(j+1)%(grid.shape[0]),k+1,6]
          grid[j,k,3]+=clonedgrid[j,k+1,7]
          grid[j,k,4]+=clonedgrid[(j-1)%(grid.shape[0]),k+1,8]
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
          grid[j,k,6]+=clonedgrid[(j-1)%(grid.shape[0]),k-1,2]
          grid[j,k,7]+=clonedgrid[j,k-1,3]
          grid[j,k,8]+=clonedgrid[(j+1)%(grid.shape[0]),k-1,4]
  return grid
  
def update(grid,relaxt,pressgradvel):
  grid=move(grid)  
  grid=relax(grid,relaxt,pressgradvel)
  return grid
  
def calc_velocity(grid):
    velocity=np.zeros((grid.shape[0],grid.shape[1],2),dtype=float)
    velocity[:,:,0]=(1/np.sum(grid,axis=2))*(grid[:,:,1]+(grid[:,:,2]+grid[:,:,8]-grid[:,:,4]-grid[:,:,6])-grid[:,:,5])    #x-direction
    velocity[:,:,1]=(1/np.sum(grid,axis=2))*(grid[:,:,3]-grid[:,:,7]+(grid[:,:,2]+grid[:,:,4]-grid[:,:,6]-grid[:,:,8]))    #y-direction
    print 'x',velocity[:,:,0]
    print 'y',velocity[:,:,1]
    plt.plot(velocity[50,:,0],np.arange(velocity.shape[1]),'r+')
    plt.show()
    return velocity
    
def add_pressure_grad(grid,pressgradvel):
    velocity=calc_velocity(grid)
    velocity[:,1:velocity.shape[1]-1,0]+=pressgradvel
    return velocity
    
def calc_eq(grid,pressgradvel):
  velocity=add_pressure_grad(grid,pressgradvel)
  Velsqr=velocity[:,:,0]**2+velocity[:,:,1]**2
  grideq=np.zeros(grid.shape,dtype=float)
  w= np.array([4./9, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36])
  e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]])
  Prod=np.dot(e.T,velocity.transpose(0,2,1))
  for i in xrange(0,9):
    grideq[:,:,i]=w[i]*np.sum(grid,axis=2)*(1+3* Prod[i] + 9/2.0* Prod[i]**2 -3.0/2 * Velsqr)
  return grideq
  
def relax(grid,relaxt,pressgradvel):
    grideq=calc_eq(grid,pressgradvel)
    grid=(1-(1/relaxt))*grid+grideq*(1/relaxt)
    return grid
    