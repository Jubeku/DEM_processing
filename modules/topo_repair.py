#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

def repairX( topoC, n_iter, threshold):
    
    # Calculate gradient of DEM
    topo_gradient = np.gradient(topoC.topo,topoC.dx,topoC.dy)
    topo_gx = topo_gradient[1]
    topo_gy = topo_gradient[0]
    
    print('Remaining iterations:', n_iter)
    n_iter -= 1 

    # Find gradients above threshold
    change = False
    topo_repair = topoC.topo.copy()
    for j in range(10, topoC.Ny-10):
        for i in range(10, topoC.Nx-10):
            if topo_gx[j,i-1] < -threshold :
                change = True
                topo_repair[j,i] =  0.9*topo_repair[j,i-1] \
                                   +0.1*topo_repair[j,i+10] 		    
            if topo_gx[j,i+1] >  threshold :
                change = True
                topo_repair[j,i] =  0.9*topo_repair[j,i+1] \
                                   +0.1*topo_repair[j,i-10] 		    
    for i in range(10, topoC.Nx-10):
        for j in range(10, topoC.Ny-10):
            if topo_gy[j-1,i] < -threshold :
                change = True
                topo_repair[j,i] =  0.9*topo_repair[j-1,i] \
                                   +0.1*topo_repair[j+10,i]
            if topo_gy[j+1,i] >  threshold :
                change = True
                topo_repair[j,i] =  0.9*topo_repair[j+1,i] \
                                   +0.1*topo_repair[j-10,i]
    topoC.topo = topo_repair

    # Recursive iteration
    if not change:
        print('DEM repaired - no more iterations needed.')
    elif n_iter > 0 :
        repairX( topoC, n_iter, threshold )	    
    else:
        print('DEM repaired.')

    return()


def repairDEM( topoC ) :
	
    """
    This function smoothes a DEM file.
    """

    print('\n Repairing ... ')

    ls = LightSource(azdeg=225, altdeg=45)
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    # Plot original DEM
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx,
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]
    ax1.set_title('Original DEM')
    ax1.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax1.set(xlabel='Easting (m)', ylabel='Northing (m)')


    # Repair DEM
    n_iter    = 20  # maximum number of recursive iterations
    threshold =  2. # threshold of gradient to repair
    repairX( topoC, n_iter, threshold )

    # Plot repaired DEM
    ax2.set_title('Repaired DEM')
    ax2.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax2.set(xlabel='Easting (m)')
    plt.show()

    return()
