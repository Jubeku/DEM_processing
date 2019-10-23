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
	
	# rockfall coordinates
	#24/04/13 rockfall
	#oe =  366245.084000
	#on = 7649885.264000
	#Stromboli crater NE
	oe =  518428.35
	on = 4293985.49
	# light for plotting topography
	ls = LightSource(azdeg=225, altdeg=45)
	cmap = plt.cm.gist_earth
	
	# read the original DEM file
	with open(pathin+fileID) as f:
	    lines     = f.readlines()
	    Nx,Ny     = map(np.int,lines[1].split())
	    east_old  = map(np.float,lines[2].split())
	    north_old = map(np.float,lines[3].split())
	    tmin,tmax = map(np.float,lines[4].split())
	    #topo      = [ map(np.float,line.split()) for line in lines[5:] ]
	    topo      = np.loadtxt(lines,skiprows=5)	
	#topo = np.array(topo)
	#topo = np.flipud(topo)
	
	# read the new coordinates
	with open(pathout+coordID) as f:
	    lines       = f.readlines()
	    east_new    = map(np.float,lines[1].split())
	    north_new   = map(np.float,lines[3].split())
            Nxout,Nyout = map(np.float,lines[5].split())
            outID       = lines[7].replace('\n', '')
            outFormat   = lines[9].replace('\n', '')
        Nxout = np.int(Nxout)
        Nyout = np.int(Nyout)

	# determine resolution of DEM file
	dx = (east_old[1]-east_old[0])/(Nx-1)
	dy = (north_old[1]-north_old[0])/(Ny-1)
	
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
	    
	    N_oe = N_oe * Nxout/Nx
	    N_on = N_on * Nyout/Ny
	    
	    Nx = Nxout
	    Ny = Nyout
		
            f, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy), cmap='gray')
	    ax1.plot(N_oe_old, N_on_old, 'r*', lw=4)	
            ax1.set_xlim([0,Nx_old])
            ax1.set_ylim([0,Ny_old])
            ax2.imshow(ls.hillshade(topo_out,vert_exag=1,dx=dxout,dy=dyout), cmap='gray')
	    ax2.plot(N_oe, N_on, 'r*', lw=4)	
            ax2.set_xlim([0,Nx])
            ax2.set_ylim([0,Ny])
            plt.show()
	else :
	    dxout = Lx/(Nxout-1)
            dyout = Ly/(Nyout-1)
            print 'Nxout: ', Nxout, ', Nyout: ', Nyout
            print 'dxout: ', dxout, ', dyout: ', dyout, '\n'

	
	if rm_mean:
	    topo_out = topo_out - np.mean(topo_out)
            f, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy), cmap='gray')
            ax1.plot(N_oe_old, N_on_old, 'r*', lw=4)
            ax1.set_xlim([0,Nx_old])
            ax1.set_ylim([0,Ny_old])
            ax2.imshow(ls.hillshade(topo_out,vert_exag=1,dx=dxout,dy=dyout), cmap='gray')
            ax2.plot(N_oe, N_on, 'r*', lw=4)
            ax2.set_xlim([0,Nx])
            ax2.set_ylim([0,Ny])
            plt.show()



	print 'Rockfall position: N_oe = ', N_oe, ', N_on = ', N_on
	
	if rm_halfhight :
	    topo_out = topo_out-(np.min(topo_out)+np.max(topo_out))/2	

        # write topo to new file
        if (outFormat=='GRD') :
            with open(pathout+outID+'.grd', 'w') as f:
                f.write('DSAA\n')
	        f.write(' '+str(Nx)+' '+str(Ny)+'\n')
	        np.savetxt(f, east_new, fmt=' %.6f', newline='')
	        f.write('\n')
	        np.savetxt(f, north_new, fmt=' %.6f', newline='')
                f.write('\n')
                np.savetxt(f, (np.min(topo_out),np.max(topo_out)), fmt=' %.1f', newline='')
                f.write('\n')
                np.savetxt(f, topo_out, fmt='%.1f', delimiter=' ')
        elif (outFormat=='SEM') :
            with open(pathout+outID+'.dat', 'w') as f:
                f.write(' ncols         '+str(Nx)+'\n')
                f.write(' nrows         '+str(Ny)+'\n')
                f.write(' xllcorner     %.6f\n' % (0.))
                f.write(' xllcorner     %.6f\n' % (0.))
                f.write(' cellsize      %.6f\n' % (dxout))
                f.write(' NODATA_value '+str(-99999)+'\n')
                np.savetxt(f, np.flipud(topo_out), fmt='   %.6f', delimiter='   ')
        elif (outFormat=='SHALTOP') :
            with open(pathout+outID+'.d', 'w') as f:
                np.savetxt(f, np.fliplr(topo_out), fmt='%.6f', delimiter='\n')
        else:
            print 'Please chose output format: GRD, SEM or SHALTOP'


	return()
