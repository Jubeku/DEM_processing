#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.interpolate import RectBivariateSpline


def taper( topoC, Nx, Ny, Nbufx, Nbufy ) :
    """
    This function adds a taper zone on each side of the DEM file.
    """

    # Original DEM dimensions and resolution
    Nx0   = topoC.Nx
    Ny0   = topoC.Ny
    dx    = topoC.dx
    dy    = topoC.dy
    # Total width and buffer width 
    Lx    = (Nx-1)*dx
    Ly    = (Ny-1)*dy
    Lbufx = Nbufx*dx
    Lbufy = Nbufy*dy
    
    # DEM with buffer
    topo_buf = np.zeros((Ny,Nx))
    topo_buf[ Nbufy:(Ny-Nbufy), Nbufx:(Nx-Nbufx) ] = topoC.topo

    # Fill buffer with border of DEM
    # edges in x-direction
    for i in range(Nbufx,(Nx-Nbufx)) :
        topo_buf[ 0:Nbufy    , i ] = topoC.topo[ 0    , i-Nbufx ]
        topo_buf[ Ny-Nbufy:Ny, i ] = topoC.topo[ Ny0-1, i-Nbufx ]
    # edges in y-direction
    for j in range(Nbufy,(Ny-Nbufy)) :
        topo_buf[ j, 0:Nbufx     ] = topoC.topo[ j-Nbufy, 0     ]
        topo_buf[ j, Nx-Nbufx:Nx ] = topoC.topo[ j-Nbufy, Nx0-1 ]
    # corners
    topo_buf[ 0:Nbufy    , 0:Nbufx     ] = topoC.topo[ 0    , 0     ]
    topo_buf[ Ny-Nbufy:Ny, 0:Nbufx     ] = topoC.topo[ Ny0-1, 0     ]
    topo_buf[ Ny-Nbufy:Ny, Nx-Nbufx:Nx ] = topoC.topo[ Ny0-1, Nx0-1 ]
    topo_buf[ 0:Nbufy    , Nx-Nbufx:Nx ] = topoC.topo[ 0    , Nx0-1 ]

    # Cosine weighting function
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

    # Weighting towards average in outward direction
    avg = np.sum(topoC.topo)/(Nx0*Ny0)
    for j in range(0,Ny) :
        y = j*dy
        wy = wtcoef(y,0.0,Lbufy,Ly-Lbufy,Ly)
        for i in range(0,Nx) :
            x = i*dx
            wx = wtcoef(x,0.0,Lbufx,Lx-Lbufx,Lx)
            topo_buf[j,i] = topo_buf[j,i]*wx*wy + avg*(1.-wx*wy)

    return(topo_buf)


def filt2d( topo, k1, k2, Nx, Ny, dx, dy, remove_avg) :

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

    topo_plot = np.zeros_like(topo)
    for j in range(0,Ny) :
        if (j <= Ny/2+1) :
            ky = j*dky
        else :
            ky = -(Ny-j)*dky
        for i in range(0,Nx) :
            if (i <= Nx/2+1) :
                kx = i*dkx
            else :
                kx = -(Nx-i)*dkx
            k = np.sqrt(kx**2+ky**2)/(2.*np.pi)
            wt = wtcoef(k,k1,k2)
            topo[j,i] = topo[j,i]*wt
            if (wt != 0) : topo_plot[j,i] = np.log10(np.abs(topo[j,i]))
    
    if remove_avg :
         topo[0,0] = 0.

    return topo


def filterDEM( topoC, k1, k2  ) :
	
    """
    This function smoothes a DEM file.
    """
    
    remove_avg = False

    # Define buffer width = 0.2 of original length on each side
    Nbufx = int(0.2*topoC.Nx)
    Nbufy = int(0.2*topoC.Ny)
    # Dimensions of DEM with buffer
    Nx    = topoC.Nx + 2*Nbufx
    Ny    = topoC.Ny + 2*Nbufy

    # Add buffer zone
    topo = taper(topoC, Nx, Ny, Nbufx, Nbufy)

    # Fourier transform to wavenumber domain
    topo = np.fft.fft2(topo)

    # Low-pass filtering with k1 and k2 corner wavenumbers
    topo = filt2d( topo, k1, k2, Nx, Ny, topoC.dx, topoC.dy, remove_avg) 

    # Inverse Fourier transform
    topo = np.real(np.fft.ifft2(topo))

    # Remove buffer zones
    topo = topo[ Nbufy:(Ny-Nbufy), Nbufx:(Nx-Nbufx) ]


    # Plotting
    ls = LightSource(azdeg=225, altdeg=45)
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    # Plot original DEM
    extent = [ topoC.E0, topoC.E0+topoC.Nx*topoC.dx,
               topoC.N0, topoC.N0+topoC.Ny*topoC.dy ]
    ax1.set_title('Original DEM')
    ax1.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax1.set(xlabel='Easting (m)', ylabel='Northing (m)')

    # Save new properties to DEM object 
    topoC.topo = topo

    # Plot cropped DEM
    ax2.set_title('Filtered DEM: '+\
            'k1='+str(k1)+'m$^{-1}$, k2='+str(k2)+'m$^{-1}$')
    ax2.imshow(ls.hillshade(topoC.topo,vert_exag=1,dx=topoC.dx,dy=topoC.dy),
               extent=extent, cmap='gray', origin='lower')
    ax2.set(xlabel='Easting (m)')
    plt.show()

    return()
