#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.interpolate import RectBivariateSpline
from modules_filter import taper, filt2d

def filter( topo, Nx, Ny, dx, dy  ) :
	
	"""
	This function smoothes a DEM file.
	"""
        
	# read filter parameters
        with open('topo_filter.dat') as f:
            lines             = f.readlines()
            Nxout,Nyout,k1,k2 = map(np.float,lines[1].split())
	    outID             = lines[3].replace('\n', '')
	    outFormat         = lines[5].replace('\n', '')
	Nxout = np.int(Nxout)
	Nyout = np.int(Nyout)

	new_grid = True	
	if (Nxout==Nx and Nyout==Ny) : new_grid = False
	
	# Define buffer width = 0.2 of original length on each side
	Nbufx = int(0.2*Nx)
	Nbufy = int(0.2*Ny)
	Lbufx = Nbufx*dx
	Lbufy = Nbufy*dy
	Nx = Nx + 2*Nbufx
	Ny = Ny + 2*Nbufy
	Lx = (Nx-1)*dx
	Ly = (Ny-1)*dy
	
	# Add buffer zone
	topo_buf = np.zeros((Ny,Nx))
	topo_buf = taper(Nx,Ny,Nbufx,Nbufy,dx,dy,topo,topo_buf)
	
	# Fourier transform to wavenumber domain
	topo_fft = np.fft.fft2(topo_buf)

	# Low-pass filtering with k1 and k2 corner wavenumbers
	topo_fft = filt2d( Nx, Ny, dx, dy, topo_fft, k1, k2, remove_avg) 
	
	# Inverse Fourier transform
	topo_filt = np.real(np.fft.ifft2(topo_fft))

        # Remove buffer zones
        topo_out = topo_filt[ Nbufy:(Ny-Nbufy), Nbufx:(Nx-Nbufx) ]
        
	# reinterpolate topo  on new grid dimensions 
	if new_grid :
	    Nx = Nx - 2*Nbufx
            Ny = Ny - 2*Nbufy	
            Lx = (Nx-1)*dx
            Ly = (Ny-1)*dy
            dxout = Lx/(Nxout-1)
            dyout = Ly/(Nyout-1)    	
            xi = np.arange(0, dx*Nx+dx, dx)
            yi = np.arange(0, dy*Ny+dy, dy)
            xo = np.arange(0, dxout*Nxout+dxout, dxout)
            yo = np.arange(0, dyout*Nyout+dyout, dyout)
            interp_spline = RectBivariateSpline(yi,xi,topo_out)
            topo_out      = interp_spline(yo,xo)
            print 'Nxout: ', Nxout, ', Nyout: ', Nyout		
            print 'dxout: ', dxout, ', dyout: ', dyout, '\n'

	
	return(topo, Nxout, Nyout, dxout, dyout)
