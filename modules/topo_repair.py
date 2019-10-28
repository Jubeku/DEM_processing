#!/usr/bin/python

import numpy as np


def repairDEM( topoC ) :
	
    """
    This function smoothes a DEM file.
    """

    # light for plotting topography
    ls = LightSource(azdeg=225, altdeg=45)
    cmap = plt.cm.gist_earth

    # read DEM file
    with open(path+fileID) as f:
        lines     = f.readlines()
        Nx,Ny     = map(np.int,lines[1].split())
        easting   = map(np.float,lines[2].split())
        northing  = map(np.float,lines[3].split())
        tmin,tmax = map(np.float,lines[4].split())
        #topo      = [ map(np.float,line.split()) for line in lines[5:] ]
        topo      = np.loadtxt(lines,skiprows=5)

    xi = np.zeros(Nx)
    yi = np.zeros(Ny) 
    for i in xrange(0,Nx) :
        xi[i] = i*dx
    for j in xrange(0,Ny) :
        yi[j] = j*dy



    topo_gradient = np.gradient(topo,0.5,0.5)
    #topo_gy = topo_gradient[0]
    #topo_gx = topo_gradient[1]
    def repairX(topo_gx,topo_gy) :
        global count
        count = count+1
        change = False
        for j in xrange(10,Ny-10):
            for i in xrange(10,Nx-10):
                if topo_gx[j,i-1] < -2 :
                    #if np.abs(topo_gx[j,i-2])-2 < 0:
                    #if topo_gx[j,i-2] < 2:
                        change = True
                        topo_repair[j,i] = 0.9*topo_repair[j,i-1]+0.1*topo_repair[j,i+10] 		    
                        #topo_repair[j,i] = 0.5*topo_repair[j,i-1]+0.5*topo_repair[j,i+1] 		    
                        #topo_repair[j,i] = (topo_repair[j+1,i]+topo_repair[j-1,i]+
                            #		topo_repair[j,i+1]+topo_repair[j,i-1])/4 		    
                        #topo_repair[j,i] = topo_repair[j,i-1] 		    
                if topo_gx[j,i+1] >  2 :
                    #if np.abs(topo_gx[j,i-2])-2 < 0:
                    #if topo_gx[j,i+2] > -2:
                        change = True
                        topo_repair[j,i] = 0.9*topo_repair[j,i+1]+0.1*topo_repair[j,i-10] 		    
                        #topo_repair[j,i] = 0.5*topo_repair[j,i-1]+0.5*topo_repair[j,i+1] 
                        #topo_repair[j,i] = (topo_repair[j+1,i]+topo_repair[j-1,i]+
                            #		topo_repair[j,i+1]+topo_repair[j,i-1])/4 		    
        for i in xrange(10,Nx-10):
            for j in xrange(10,Ny-10):
                if topo_gy[j-1,i] < -2 :
                    #if np.abs(topo_gx[j,i-2])-2 < 0:
                    #if topo_gx[j,i-2] < 2:
                        change = True	
                        topo_repair[j,i] = 0.9*topo_repair[j-1,i]+0.1*topo_repair[j+10,i]
                        #topo_repair[j,i] = (topo_repair[j+1,i]+topo_repair[j-1,i]+
                            #		topo_repair[j,i+1]+topo_repair[j,i-1])/4 		    
                        #topo_repair[j,i] = topo_repair[j,i-1] 		    
                if topo_gy[j+1,i] >  2 :
                    #if np.abs(topo_gx[j,i-2])-2 < 0:
                    #if topo_gx[j,i+2] > -2:
                        change = True
                        topo_repair[j,i] = 0.9*topo_repair[j+1,i]+0.1*topo_repair[j-10,i]
                        #topo_repair[j,i] = 0.5*topo_repair[j,i-1]+0.5*topo_repair[j,i+1] 
                        #topo_repair[j,i] = (topo_repair[j+1,i]+topo_repair[j-1,i]+
                            #		topo_repair[j,i+1]+topo_repair[j,i-1])/4 		    
        print(change)
        if change == True and count < 20 :
            topo_gradient_repair = np.gradient(topo_repair,0.5,0.5)
            repairX(topo_gradient_repair[1],topo_gradient_repair[0])	    
        return()
            
    topo_repair = topo.copy()
    global count 
    count = 0
    repairX(topo_gradient[1],topo_gradient[0])	
    topo_gradient_repair = np.gradient(topo_repair,0.5,0.5)



    return()
