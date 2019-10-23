#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.interpolate import RectBivariateSpline
from modules_filter import taper, filt2d

def filter( path, fileID, parameterID  ) :
	
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
        #topo = np.array(topo)
        topo = np.flipud(topo)

	topo[topo < 0.] = 0.0
        
	# read filter parameters
        with open(path+parameterID) as f:
            lines             = f.readlines()
            Nxout,Nyout,k1,k2 = map(np.float,lines[1].split())
	    outID             = lines[3].replace('\n', '')
	    outFormat         = lines[5].replace('\n', '')
	Nxout = np.int(Nxout)
	Nyout = np.int(Nyout)
	new_grid = True	
	if (Nxout==Nx and Nyout==Ny) : new_grid = False
		
        # determine resolution of DEM file
        dx = (easting[1]-easting[0])/(Nx-1)
        dy = (northing[1]-northing[0])/(Ny-1)
	print 'Nx: ', Nx, ', Ny: ', Ny 
	print 'dx: ', dx, ', dy: ', dy, '\n'
	
	# define buffer width = 0.2 of original length on each side
	Nbufx = int(0.2*Nx)
	Nbufy = int(0.2*Ny)
	Lbufx = Nbufx*dx
	Lbufy = Nbufy*dy
	Nx = Nx + 2*Nbufx
	Ny = Ny + 2*Nbufy
	Lx = (Nx-1)*dx
	Ly = (Ny-1)*dy
	
	print np.shape(topo)
	
	# add buffer zone
	topo_buf = np.zeros((Ny,Nx))
	topo_buf = taper(Nx,Ny,Nbufx,Nbufy,dx,dy,topo,topo_buf)

	xi = np.zeros(Nx)
	yi = np.zeros(Ny) 
	for i in xrange(0,Nx) :
	    xi[i] = i*dx
	for j in xrange(0,Ny) :
	    yi[j] = j*dy
        
	f, (ax1, ax2) = plt.subplots(1, 2)
        #ax1.imshow(topo, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
        #ax2.imshow(topo_buf, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
        ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy), extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	ax2.imshow(ls.hillshade(topo_buf,vert_exag=1,dx=dx,dy=dy), extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	plt.show()
	
	"""
	# nearest number of power of 2 for Nx and Ny
	Nxp = 2**(int(np.log(float(Nx))/np.log(2.0))+1)
	Nyp = 2**(int(np.log(float(Ny))/np.log(2.0))+1)
	print Nxp, Nyp
	print('Filtering...')


	# reinterpolate topography grid by spline on the FFT points
	dxp=Lx/(Nxp-1)
	dyp=Ly/(Nyp-1)
 	
	topo_fft = np.zeros((Nxp,Nyp))
	xi = np.zeros(Nx)
	yi = np.zeros(Ny) 
	xp = np.zeros(Nxp)
	yp = np.zeros(Nyp)

	for i in xrange(0,Nx) :
	    xi[i] = i*dx
	for j in xrange(0,Ny) :
	    yi[j] = j*dy
	for i in xrange(0,Nxp) :
	    xp[i] = i*dxp
	for j in xrange(0,Nyp) :
	    yp[j] = j*dyp

	print Nx, Ny, Nxp, Nyp

	
	xxi, yyi = np.meshgrid(xi, yi)
	xxp, yyp = np.meshgrid(xp, yp)

	#topo_rbf_func = interp.Rbf(xxi,yyi,topo_buf,smooth=0)
	#topo_rbf      = topo_rbf_func(xxp,yyp)
	
	interp_spline = RectBivariateSpline(xi,yi,topo_buf)
	topo_rbf      = interp_spline(xp,yp)

	f, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(topo_buf, aspect='auto', vmin=tmin, vmax=tmax)
        ax2.imshow(topo_rbf, aspect='auto', vmin=tmin, vmax=tmax)
        plt.show()	

	"""

	# Fourier transform to wavenumber domain
	#topo_fft = np.fft.fft2(topo_rbf)	
	topo_fft = np.fft.fft2(topo_buf)

	# low-pass filtering with k1 and k2 corner wavenumbers
	remove_avg = True 
	#topo_fft = filt2d( Nxp, Nyp, dxp, dyp, topo_fft, k1, k2, remove_avg) 
	topo_fft = filt2d( Nx, Ny, dx, dy, topo_fft, k1, k2, remove_avg) 
	
	# inverse Fourier transform
	topo_filt = np.real(np.fft.ifft2(topo_fft))

        f, (ax1, ax2) = plt.subplots(1, 2)
        #ax1.imshow(topo_buf, aspect='auto', extent=[xi[0],xi[-1],yi[0],yi[-1]], vmin=tmin, vmax=tmax)
        #ax2.imshow(topo_filt, aspect='auto', extent=[xi[0],xi[-1],yi[0],yi[-1]], vmin=tmin, vmax=tmax)
	ax1.imshow(ls.hillshade(topo_buf,vert_exag=1,dx=dx,dy=dy), 
		   extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
        ax2.imshow(ls.hillshade(topo_filt,vert_exag=1,dx=dx,dy=dy), 
		   extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	ax1.set_title('original topography')
	ax2.set_title('low pass: k1=0.003, k2=0.005')
	ax1.set_ylabel('North (m)')
	ax1.set_xlabel('East (m)')
	ax2.set_xlabel('East (m)')
	plt.show()


        # remove buffer zones
        topo_out = topo_filt[ Nbufy:(Ny-Nbufy), Nbufx:(Nx-Nbufx) ]
        
	# reinterpolate topo  on new grid dimensions 
	if new_grid :
	    Nx = Nx - 2*Nbufx
            Ny = Ny - 2*Nbufy	
            Lx = (Nx-1)*dx
            Ly = (Ny-1)*dy
            dxout = Lx/(Nxout-1)
            dyout = Ly/(Nyout-1)    	
            xi = np.zeros(Nx)
            yi = np.zeros(Ny)
	    xo = np.zeros(Nxout)
            yo = np.zeros(Nyout)
            for i in xrange(0,Nx) :
                xi[i] = i*dx
            for j in xrange(0,Ny) :
                yi[j] = j*dy
            for i in xrange(0,Nxout) :
                xo[i] = i*dxout
            for j in xrange(0,Nyout) :
                yo[j] = j*dyout
            interp_spline = RectBivariateSpline(yi,xi,topo_out)
            topo_out      = interp_spline(yo,xo)
            print 'Nxout: ', Nxout, ', Nyout: ', Nyout		
            print 'dxout: ', dxout, ', dyout: ', dyout, '\n'

        f, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(topo, aspect='auto', vmin=tmin, vmax=tmax)
        ax2.imshow(topo_out, aspect='auto', vmin=tmin, vmax=tmax)
        plt.show()      


	
	return()
