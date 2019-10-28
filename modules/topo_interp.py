#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline

def interpDEM( topoC, dxout, dyout ) :

    """
    This function takes a DEM file in fileID (format .grd) and cuts it to a 
    smaller	region with coordiates specified in parameter file coordID.
    """
   
    print('\n Interpolating ... ')

    # Dimensions of DEM
    Lx = (topoC.Nx-1)*topoC.dx
    Ly = (topoC.Ny-1)*topoC.dy
    # New number of grid points
    Nxout = np.int(Lx/dxout) + 1 
    Nyout = np.int(Ly/dyout) + 1

    # Define coordinate arrays
    xi = np.arange(0., Lx+topoC.dx, topoC.dx)
    yi = np.arange(0., Ly+topoC.dy, topoC.dy)
    xo = np.arange(0., Lx+dxout, dxout)
    yo = np.arange(0., Ly+dyout, dyout)
    
    # Interpolation on new grid coordinates
    f_interp = RectBivariateSpline(yi, xi, topoC.topo)
    topo     = f_interp(yo,xo)
    
    
    print('Nx: ', topoC.Nx, ', Ny: ', topoC.Ny)
    print('dx: ', topoC.dx, ', dy: ', topoC.dy)
    print('Nxout: ', Nxout, ', Nyout: ', Nyout)
    print('dxout: ', dxout, ', dyout: ', dyout)
    
    # Plotting
    ls = LightSource(azdeg=225, altdeg=45)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx,
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]

    # Plot original DEM
    ax1.set_title('Original')
    ax1.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax1.set(xlabel='Easting (m)', ylabel='Northing (m)')

    # Save new properties to DEM object
    topoC.topo = topo
    topoC.Nx   = Nxout
    topoC.Ny   = Nyout

    # Plot cropped DEM
    ax2.set_title('Interpolated')
    ax2.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax2.set(xlabel='Easting (m)')
    plt.show()

    return()
