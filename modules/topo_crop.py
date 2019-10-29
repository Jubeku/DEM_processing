#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

def cropDEM( topoC, E0out, E1out, N0out, N1out ) :

    """
    This function crops a DEM object to new coordinates. 
    
    Input paramters:
    topoC  :  DEM object of class Topo
    E0out  :  Left easting coordinate
    E1out  :  Right easting coordinate
    N0out  :  Bottom northing coordinate
    N1out  :  Top northing coordinate
 
    """
    
    print('\n Cropping ... ')

    # Determine indeces which correspond to new coordinates, be aware 
    # of the grd format convention: ordered from left to right (Easting) 
    # and bottom to top (Northing)
    idx_e0 = np.int((E0out-topoC.E0)/topoC.dx)
    idx_e1 = np.int((E1out-topoC.E0)/topoC.dx)+1
    idx_n0 = np.int((N0out-topoC.N0)/topoC.dy)
    idx_n1 = np.int((N1out-topoC.N0)/topoC.dy)+1
    
    # New grid dimensions
    Nx = idx_e1-idx_e0
    Ny = idx_n1-idx_n0
        
    # Crop DEM to new coordinates
    topo = np.zeros((Ny,Nx)) 
    topo = topoC.topo[idx_n0:idx_n1,idx_e0:idx_e1]
    
    # Plotting
    ls = LightSource(azdeg=225, altdeg=45) 
    fig, (ax1, ax2) = plt.subplots(1, 2)
   
    # Plot original DEM
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx,
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]
    ax1.set_title('Original DEM')
    ax1.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy), 
               extent=extent, cmap='gray', origin='lower')
    ax1.plot([E0out, E1out], [N0out, N0out], 'k-', lw=2)	
    ax1.plot([E1out, E1out], [N0out, N1out], 'k-', lw=2)	
    ax1.plot([E1out, E0out], [N1out, N1out], 'k-', lw=2)	
    ax1.plot([E0out, E0out], [N1out, N0out], 'k-', lw=2)	
    ax1.set(xlabel='Easting (m)', ylabel='Northing (m)') 

    # Save new properties to DEM object 
    topoC.topo = topo
    topoC.E0   = E0out
    topoC.N0   = N0out
    topoC.Nx   = Nx
    topoC.Ny   = Ny
    
    # Plot cropped DEM
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx,
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]
    ax2.set_title('Cropped DEM')
    ax2.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy), 
               extent=extent, cmap='gray', origin='lower')
    ax2.set(xlabel='Easting (m)') 
    plt.show()

    return()
