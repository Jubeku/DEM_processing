#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

def taper( Nx, Ny, Nbufx, Nbufy, dx, dy, topo, topo_buf ) :
    """
    This function adds a taper zone on each side of the DEM file.
    """
    
    Nx0 = Nx - 2*Nbufx
    Ny0 = Ny - 2*Nbufy
    Lbufx = Nbufx*dx
    Lbufy = Nbufy*dy
    Lx = (Nx-1)*dx
    Ly = (Ny-1)*dy

    avg = np.sum(topo)/(Nx0*Ny0)
    topo_buf[ Nbufy:(Ny-Nbufy), Nbufx:(Nx-Nbufx) ] = topo
    
    # edges in x-direction
    for i in xrange(Nbufx,(Nx-Nbufx)) :
        topo_buf[ 0:Nbufy    , i ] = topo[ 0    , i-Nbufx ]
        topo_buf[ Ny-Nbufy:Ny, i ] = topo[ Ny0-1, i-Nbufx ]
    # edges in y-direction
    for j in xrange(Nbufy,(Ny-Nbufy)) :
        topo_buf[ j, 0:Nbufx     ] = topo[ j-Nbufy, 0     ]
        topo_buf[ j, Nx-Nbufx:Nx ] = topo[ j-Nbufy, Nx0-1 ]
    # corners
    topo_buf[ 0:Nbufy    , 0:Nbufx     ] = topo[ 0    , 0     ]
    topo_buf[ Ny-Nbufy:Ny, 0:Nbufx     ] = topo[ Ny0-1, 0     ]
    topo_buf[ Ny-Nbufy:Ny, Nx-Nbufx:Nx ] = topo[ Ny0-1, Nx0-1 ]
    topo_buf[ 0:Nbufy    , Nx-Nbufx:Nx ] = topo[ 0    , Nx0-1 ]
    

    def wtcoef(p, p0, p1, p2, p3) :
        if p1 <= p <= p2 :
            wt = 1.0
        elif p0 <= p < p1 :
            wt = 0.5 * (1. + np.cos(np.pi*(p-p1)/(p1-p0)) )
        elif p2 < p <= p3 :
            wt = 0.5 * (1. + np.cos(np.pi*(p-p2)/(p3-p2)) )
        else :
            print('Please check buffer coordinates!')
        return wt 
   
    # weighting towards average in cosine shape
    for j in xrange(0,Ny) :
        y = j*dy
        wy = wtcoef(y,0.0,Lbufy,Ly-Lbufy,Ly) 
        for i in xrange(0,Nx) :
            x = i*dx
            wx = wtcoef(x,0.0,Lbufx,Lx-Lbufx,Lx)
            topo_buf[j,i] = topo_buf[j,i]*wx*wy + avg*(1.-wx*wy)

    return topo_buf


def filt2d( Nx, Ny, dx, dy, topo, k1, k2, remove_avg) :
	
    """ 
    This function is a low-pass on topo in the wavemumber domain. k1 and k2 are 
    the corner filter coefficients with cosinal transition in between.
    """
    topo_nofilt = topo.copy()

    def wtcoef(k, k1, k2) :
        if k1 > k2 :
            print('check filter coefficients: k1 > k2 :(')
        elif k >= k2 :
            wt = 0.0
        elif k < k1 :
            wt = 1.0
        else :
            wt = 0.5 * (1.0 + np.cos(np.pi*(k-k1)/(k2-k1)) )
        return wt

    dkx = 2.*np.pi/np.float(Nx-1)/dx
    dky = 2.*np.pi/np.float(Ny-1)/dy
    print 'dkx: ', dkx, ', dky: ', dky 
    print 'kxmax: ', dkx*(Nx/2+1), ', kymax: ', dky*(Ny/2+1), '\n' 
    
    topo_plot = np.zeros_like(topo)	
    for j in xrange(0,Ny) :
        if (j <= Ny/2+1) :
            ky = j*dky
        else :
            ky = -(Ny-j)*dky
        for i in xrange(0,Nx) :
            if (i <= Nx/2+1) :
                kx = i*dkx
            else :
                kx = -(Nx-i)*dkx
            k = np.sqrt(kx**2+ky**2)/(2.*np.pi)
            wt = wtcoef(k,k1,k2)
            topo[j,i] = topo[j,i]*wt
            if (wt != 0) : topo_plot[j,i] = np.log10(np.abs(topo[j,i]))
    
    #plotting:
    kxi_log = np.logspace(np.log10(dkx),np.log10(dkx*(Nx/2+1)),Nx/2)
    kyi_log = np.logspace(np.log10(dky),np.log10(dky*(Ny/2+1)),Ny/2)
    X, Y = np.meshgrid(kxi_log,kyi_log)
    
    f, (ax1, ax2) = plt.subplots(1, 2)
    ax1.pcolor(X,Y,np.log10(np.abs(topo_nofilt[1:Ny/2+1,1:Nx/2+1])))#, vmin=tmin, vmax=tmax)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([kxi_log.min(), kxi_log.max()])
    ax1.set_ylim([kyi_log.min(), kyi_log.max()])
    ax1.grid(True)
    ax2.pcolor(X,Y,topo_plot[1:Ny/2+1,1:Nx/2+1])#, vmin=tmin, vmax=tmax)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim([kxi_log.min(), kxi_log.max()])
    ax2.set_ylim([kyi_log.min(), kyi_log.max()])
    ax2.grid(True)
    ax1.set_title('original topography')
    ax2.set_title('low pass: k1=0.003, k2=0.005')
    ax1.set_ylabel('ky (1/m)')
    ax1.set_xlabel('kx (1/m)')
    ax2.set_xlabel('kx (1/m)')
    plt.show()
    
    if remove_avg :
         print topo[0,0]
         topo[0,0] = 0.

    return topo


	
