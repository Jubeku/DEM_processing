#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline

def cut( pathin, fileID, pathout, coordID, rm_halfhight, rm_mean ) :

	"""
	This function takes a DEM file in fileID (format .grd) and cuts it to a 
	smaller	region with coordiates specified in parameter file coordID.
	"""
	
	# read the new coordinates
	with open(pathout+coordID) as f:
	    lines       = f.readlines()
	    east_new    = map(np.float,lines[1].split())
	    north_new   = map(np.float,lines[3].split())
            Nxout,Nyout = map(np.float,lines[5].split())
        Nxout = np.int(Nxout)
        Nyout = np.int(Nyout)

	# determine indeces which correspond to new coordinates, be aware 
	# of the grd format convention: ordered from left to right (Easting) 
	# and bottom to top (Northing)
	idx_e0 = np.int((east_new[0]-east_old[0])/dx)
	idx_e1 = np.int((east_new[1]-east_old[0])/dx)+1
	idx_n0 = np.int((north_new[0]-north_old[0])/dy)
	idx_n1 = np.int((north_new[1]-north_old[0])/dy)+1

	Nx_old = Nx	
	Ny_old = Ny	
	Nx = idx_e1-idx_e0
	Ny = idx_n1-idx_n0
        Lx = (Nx-1)*dx
        Ly = (Ny-1)*dy
        new_grid = True
        if (Nxout==Nx and Nyout==Ny) : new_grid = False
	print 'Lx (new): ', Lx, ', Ly (new): ', Ly
	print 'Nx (new): ', Nx, ', Ny (new): ', Ny
	print 'Interpolation: ', new_grid
	
	# cut file to new coordinates
	print np.shape(topo)
	topo_out = np.zeros((Ny,Nx)) 
	topo_out = topo[idx_n0:idx_n1,idx_e0:idx_e1]

	tmin = 2200
	N_oe_old = np.int((oe-east_old[0])/dx)
	N_on_old = np.int((on-north_old[0])/dy)
	N_oe = N_oe_old-idx_e0
	N_on = N_on_old-idx_n0
	f, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy), cmap='gray')
	ax1.plot([idx_e0, idx_e1], [idx_n0, idx_n0], 'k-', lw=2)	
	ax1.plot([idx_e1, idx_e1], [idx_n0, idx_n1], 'k-', lw=2)	
	ax1.plot([idx_e1, idx_e0], [idx_n1, idx_n1], 'k-', lw=2)	
	ax1.plot([idx_e0, idx_e0], [idx_n1, idx_n0], 'k-', lw=2)	
	ax1.plot(N_oe_old, N_on_old, 'r*', lw=4)	
        ax1.set_xlim([0,Nx_old])
        ax1.set_ylim([0,Ny_old])
	ax2.imshow(ls.hillshade(topo_out,vert_exag=1,dx=dx,dy=dy), cmap='gray')
	ax2.plot(N_oe, N_on, 'r*', lw=4)	
        ax2.set_xlim([0,Nx])
        ax2.set_ylim([0,Ny])
	plt.show()


	X = np.arange(0.0,Nx*dx,dx)
	Y = np.arange(0.0,Ny*dy,dy)
	X, Y = np.meshgrid(X, Y)
	Z = topo_out.reshape(X.shape)
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	#rgb = ls.shade(Z, cmap='gray', vert_exag=1, blend_mode='soft')
	ax.plot_surface(X,Y,Z,linewidth=0, antialiased=False, shade=True) 


	
        # reinterpolate topo  on new grid dimensions 
        if new_grid :
            dxout = Lx/(Nxout-1)
            dyout = Ly/(Nyout-1)
            xi = np.arange(0, dx*Nx+dx, dx)
            yi = np.arange(0, dy*Ny+dy, dy)
            xo = np.arange(0, dxout*Nxout+dxout, dxout)
            yo = np.arange(0, dyout*Nyout+dyout, dyout)
            interp_spline = RectBivariateSpline(yi,xi,topo_out)
            topo_out      = interp_spline(yo,xo)
	    N_oe = N_oe * Nxout/Nx
	    N_on = N_on * Nyout/Ny
	    Nx = Nxout
	    Ny = Nyout
	else :
	    dxout = Lx/(Nxout-1)
            dyout = Ly/(Nyout-1)
        print('Nxout: ', Nxout, ', Nyout: ', Nyout)
        print('dxout: ', dxout, ', dyout: ', dyout)


	return()
