#!/usr/bin/python

import numpy as np
import scipy.linalg
import obspy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from scipy.interpolate import RectBivariateSpline


def repair( path, fileID  ) :
	
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

        # determine resolution of DEM file
        dx = (easting[1]-easting[0])/(Nx-1)
        dy = (northing[1]-northing[0])/(Ny-1)
	print 'Nx: ', Nx, ', Ny: ', Ny 
	print 'dx: ', dx, ', dy: ', dy, '\n'
	
	xi = np.zeros(Nx)
	yi = np.zeros(Ny) 
	for i in xrange(0,Nx) :
	    xi[i] = i*dx
	for j in xrange(0,Ny) :
	    yi[j] = j*dy
      
	np.savetxt(path+'lala.txt',topo[481:501,100:120], fmt='%.1f', delimiter=' ')
	
	"""
	X, Y = np.meshgrid(xi,yi)
	XX = X.flatten()
	YY = X.flatten()
	
	data = np.reshape(topo,XX.shape)
	
	A = np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY*2]
	C,_,_,_ = scipy.linalg.lstsq(A, data)

	Z = np.dot(A, C).reshape(X.shape)

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
	plt.xlabel('X')
	plt.ylabel('Y')
	ax.set_zlabel('Z')
	ax.axis('equal')
	ax.axis('tight')
	plt.show()
	"""

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
	    print change
	    if change == True and count < 20 :
	        topo_gradient_repair = np.gradient(topo_repair,0.5,0.5)
	        """
		f, (ax1, ax2) = plt.subplots(1, 2)
	        im1 = ax1.imshow(topo_gradient[1],vmin=-5,vmax=5)
	        f.colorbar(im1, ax=ax1)
	        im2 = ax2.imshow(topo_gradient_repair[1],vmin=-5,vmax=5)
	        f.colorbar(im2, ax=ax2)
	        plt.show()

                f, (ax1, ax2) = plt.subplots(1, 2)
                im1 = ax1.imshow(topo_gradient[0],vmin=-5,vmax=5)
                f.colorbar(im1, ax=ax1)
                im2 = ax2.imshow(topo_gradient_repair[0],vmin=-5,vmax=5)
                f.colorbar(im2, ax=ax2)
                plt.show()
		"""
		"""
		f, (ax1, ax2) = plt.subplots(1, 2)
	        #ax1.imshow(topo, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
	        #ax2.imshow(topo_buf, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
	        ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy),
        	                        extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	        ax2.imshow(ls.hillshade(topo_repair,vert_exag=1,dx=dx,dy=dy),
                	                extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	        plt.show()
		"""
		repairX(topo_gradient_repair[1],topo_gradient_repair[0])	    
	    return()
		
	topo_repair = topo.copy()
	global count 
	count = 0
        repairX(topo_gradient[1],topo_gradient[0])	
	topo_gradient_repair = np.gradient(topo_repair,0.5,0.5)

		    #if (topo[j,i-1]-topo[j,i])>2 and (topo[j,i+2]-topo[j,i])>1 :
			#topo_repair[j,i] = (topo[j,i-1]+topo[j,i+2])/2
	 	    #if (topo[j,i-1]-topo[j,i])>2
		    #if (-topo_repair[j,i-1]+topo_repair[j,i])>2 :
			#topo_repair[j,i] = topo_repair[j,i-1]


	f, (ax1, ax2) = plt.subplots(1, 2)
        #ax1.imshow(topo, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
        #ax2.imshow(topo_buf, extent=[xi[0],xi[-1],yi[0],yi[-1]], aspect='auto')#, vmin=tmin, vmax=tmax)
        ax1.imshow(ls.hillshade(topo,vert_exag=1,dx=dx,dy=dy), 
				extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	ax2.imshow(ls.hillshade(topo_repair,vert_exag=1,dx=dx,dy=dy), 
				extent=[xi[0],xi[-1],yi[0],yi[-1]], cmap='gray')
	plt.show()
	
	f, (ax1, ax2) = plt.subplots(1, 2)
	im1 = ax1.imshow(topo_gradient[1],vmin=-10,vmax=10)
	f.colorbar(im1, ax=ax1)
	im2 = ax2.imshow(topo_gradient_repair[1],vmin=-10,vmax=10)
	f.colorbar(im2, ax=ax2)
	plt.show()
	
	f, (ax1, ax2) = plt.subplots(1, 2)
	im1 = ax1.imshow(topo_gradient[0],vmin=-10,vmax=10)
	f.colorbar(im1, ax=ax1)
	im2 = ax2.imshow(topo_gradient_repair[0],vmin=-10,vmax=10)
	f.colorbar(im2, ax=ax2)
	plt.show()
		
	
        # write topo to new file
        with open(path+fileID+'_repaired.grd', 'w') as f:
            f.write('DSAA\n')
            f.write(' '+str(Nx)+' '+str(Ny)+'\n')
            np.savetxt(f, easting, fmt=' %.6f', newline='')
            f.write('\n')
            np.savetxt(f, northing, fmt=' %.6f', newline='')
            f.write('\n')
            np.savetxt(f, (np.min(topo_repair),np.max(topo_repair)), fmt=' %.1f', newline='')
            f.write('\n')
            np.savetxt(f, topo_repair, fmt='%.1f', delimiter=' ')
	
	
	return()
