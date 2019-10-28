#!/usr/bin/env python
"""
This script can be used to repair, filter, and crop 2d DEM files. 
"""

import numpy as np
import matplotlib.pyplot as plt
from modules.topo_filter import filterDEM
from modules.topo_repair import repairDEM
from modules.topo_crop import cropDEM


def main():
    
    ### INPUT
    # read input 
    with open('input.txt') as f:
        lines       = f.readlines()
        repair_bool = lines[1].replace('\n', '')
        filt_bool   = lines[3].replace('\n', '')
        cut_bool    = lines[5].replace('\n', '')
        fileID      = lines[7].replace('\n', '')
        outID       = lines[9].replace('\n', '')
    
    # read DEM file
    with open(fileID+'.grd') as f:
        lines     = f.readlines()
        Nx,Ny     = map(np.int,lines[1].split())
        easting   = map(np.float,lines[2].split())
        northing  = map(np.float,lines[3].split())
        tmin,tmax = map(np.float,lines[4].split())
        topo      = np.loadtxt(lines,skiprows=5)
    
    topo = np.flipud(topo)
    # determine resolution of DEM file
    dx = (easting[1]-easting[0])/(Nx-1)
    dy = (northing[1]-northing[0])/(Ny-1)
    xi = np.arange(0, dx*Nx+dx, dx)
    yi = np.arange(0, dy*Ny+dy, dy)
    print('Nx: ', Nx, ', Ny: ', Ny)
    print('dx: ', dx, ', dy: ', dy)

    topo_raw = topo.copy() # save copy for comparison 
    
    
    ### PROCESSING
    if repair_bool:
        topo = repairDEM(topo)
    if filter_bool:
        topo = filterDEM(topo)
    if cut_bool:
        topo = cropDEM(topo)
    
    ### PLOTTING
    ls = LightSource(azdeg=225, altdeg=45)
    #cmap = plt.cm.gist_earth
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.imshow(ls.hillshade(topo_raw,vert_exag=1,dx=dx,dy=dy),
               extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
    ax2.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy),
               extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
    ax1.set_title('original topography')
    ax2.set_title('low pass: k1=0.003, k2=0.005')
    ax1.set_ylabel('North (m)')
    ax1.set_xlabel('East (m)')
    ax2.set_xlabel('East (m)')
    plt.show()
    

    ### WRITING
    with open(outID+'.grd', 'w') as f:
        f.write('DSAA\n')
        f.write(' '+str(Nxout)+' '+str(Nyout)+'\n')
        np.savetxt(f, easting, fmt=' %.6f', newline='')
        f.write('\n')
        np.savetxt(f, northing, fmt=' %.6f', newline='')
        f.write('\n')
        np.savetxt(f, (np.min(topo_out),np.max(topo_out)), fmt=' %.1f', newline='')
        f.write('\n')
        np.savetxt(f, topo_out, fmt='%.3f', delimiter=' ')


if __name__ == "__main__":
    main()
