#!/usr/bin/env python
"""
This script can be used to repair, filter, and crop 2d DEM files. 
"""

import numpy as np
import matplotlib.pyplot as plt
from modules.classTopo import Topo
from modules.m_PLOTS import plotDEM


def main():
   

    ### INPUT
    # Read input parameter 
    with open('input.txt') as f:
        lines        = f.readlines()
        repair_bool  = lines[1].replace('\n', '')
        k1, k2i      = map(np.float,lines[3].split())
        E0out, E1out = map(np.float,lines[5].split())
        N0out, N1out = map(np.float,lines[7].split())
        dxout, dyout = map(np.float,lines[9].split())
        fileID       = lines[11].replace('\n', '')
        outID        = lines[13].replace('\n', '')
    
    # Read DEM file
    with open(fileID+'.grd') as f:
        lines     = f.readlines()
        Nx, Ny    = map(np.int,lines[1].split())
        E0, E1    = map(np.float,lines[2].split())
        N0, N1    = map(np.float,lines[3].split())
        tmin,tmax = map(np.float,lines[4].split())
        topo      = np.loadtxt(lines,skiprows=5)
   
    # Determine resolution of DEM file
    dx = (E1-E0)/(Nx-1)
    dy = (N1-N0)/(Ny-1)
    xi = np.arange(0, dx*Nx+dx, dx)
    yi = np.arange(0, dy*Ny+dy, dy)
    print('Nx: ', Nx, ', Ny: ', Ny)
    print('dx: ', dx, ', dy: ', dy)

    # Creat object with Topo class 
    topo = np.flipud(topo)
    topoC = Topo(topo, E0, N0, dx, dy, Nx, Ny)

    ### PROCESSING
    # Repairing
    if repair_bool == 'True':
        topoC.repair()
    # Filtering
    if k1 == 0.:
        print('No filtering')
    else:
        topoC.filter( k1, k2 )
    # Cropping
    if ( E0out == E0 and E1out == E1 and N0out == N0 and N1out == N1 ): 
        print('No cropping')
    else:
        topoC.crop( E0out, E1out, N0out, N1out )
    # Interpolating
    if ( dxout == dx and dyout == dy ): 
        print('No interpolation')
    else:
        topoC.interpolate( dxout, dyout )

    ### PLOTTING
    topoC.plot()

    ### WRITING
    with open(outID+'.grd', 'w') as f:
        f.write('DSAA\n')
        f.write(' '+str(topoC.Nx)+' '+str(topoC.Ny)+'\n')
        f.write(' '+str(E0out)+' '+str(E1out)+'\n')
        f.write(' '+str(N0out)+' '+str(N1out)+'\n')
        np.savetxt(f,(np.min(topoC.topo),np.max(topoC.topo)),
                fmt=' %.1f',newline='')
        f.write('\n')
        np.savetxt(f, np.flipud(topoC.topo), fmt='%.3f', delimiter=' ')


if __name__ == "__main__":
    main()
