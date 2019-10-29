#!/usr/bin/env python

from modules.topo_repair import repairDEM
from modules.topo_filter import filterDEM
from modules.topo_crop import cropDEM
from modules.topo_interp import interpDEM
from modules.m_PLOTS import plotDEM

class Topo:
    
    """
    Objects of this class represent a Digital Elevation Model (DEM).
    

    Functions and parameters:

    __init__: Object initialization.
    -------------------------------------------------------
    topo:  Array of size (Nx, Ny) of topographic elevation
    E0  :  Easting value defining lower left corner of DEM
    N0  :  Northing value defining lower left corner of DEM
    dx  :  Resolution step size in x-direction (East)
    dy  :  Resolution step size in y-direction (North)
    Nx  :  Grid dimension in x-direction
    Ny  :  Grid dimension in y-direction

    repair:  Repairing function

    filter:  Filtering function
    -------------------------------------------------------
    k1 : Corner wavenumber
    k2 : Maximum wavenumber
    
    crop:  Cropping function
    -------------------------------------------------------
    E0out  :  Left easting coordinate
    E1out  :  Right easting coordinate
    N0out  :  Bottom northing coordinate
    N1out  :  Top northing coordinate
    
    interpolate:  Interplating function
    -------------------------------------------------------
    dxout  :  Nre resolution step size in x-direction (East)
    dyout  :  New resolution step size in y-direction (North)
    
    plot:  Interplating function
    -------------------------------------------------------
    title  :  Figure title (string)
    """

    def __init__(self, topo, E0, N0, dx, dy, Nx, Ny):
        self.topo = topo
        self.E0 = E0
        self.N0 = N0
        self.dx = dx
        self.dy = dy
        self.Nx = Nx
        self.Ny = Ny

    def repair(self):
        repairDEM(self)
    
    def filter( self, k1, k2 ):
        filterDEM( self, k1, k2 )

    def crop( self, E0out, E1out, N0out, N1out ):
        cropDEM(self, E0out, E1out, N0out, N1out )
    
    def interpolate( self, dxout, dyout ):
        interpDEM(self, dxout, dyout )

    def plot(self, title):
        plotDEM(self, title)
