#!/usr/bin/env python
"""
This script demonstrates how to use matplotlib for 2d or 3d colormaps
"""

import numpy as np
import matplotlib.pyplot as plt
from topo_filter import filterDEM
from topo_repair import repairDEM
from topo_cut import cutDEM


def main():

    # read input 
    with open('input.txt') as f:
        lines       = f.readlines()
        repair_bool = lines[1].replace('\n', '')
        filt_bool   = lines[3].replace('\n', '')
        cut_bool    = lines[5].replace('\n', '')
        fileID      = lines[7].replace('\n', '')
        outID       = lines[9].replace('\n', '')
    
    # read DEM file
    with open(fileID) as f:
        lines     = f.readlines()
        Nx,Ny     = map(np.int,lines[1].split())
        easting   = map(np.float,lines[2].split())
        northing  = map(np.float,lines[3].split())
        tmin,tmax = map(np.float,lines[4].split())
        topo      = np.loadtxt(lines,skiprows=5)
    topo = np.flipud(topo)

    if repair_bool:

        repairDEM()

    if filter_bool:

        filterDEM()

    if cut_bool:

        cutDEM()
    

    
    # write topo to new file
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
