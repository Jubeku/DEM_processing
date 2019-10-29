### This defines a class 

from modules.topo_repair import repairDEM
from modules.topo_filter import filterDEM
from modules.topo_crop import cropDEM
from modules.topo_interp import interpDEM
from modules.m_PLOTS import plotDEM

class Topo:
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
