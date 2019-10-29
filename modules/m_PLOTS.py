#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

def plotDEM( topoC, title ) :

    ls = LightSource(azdeg=225, altdeg=45) 
    
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx, 
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_title( title)
    ax.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
              extent=extent, cmap='gray', origin='lower')
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    
    plt.show()
    
    return()

