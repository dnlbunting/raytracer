#!/usr/bin/env python
from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from .optical_elements import *
from .sources import Source, SingleRay, CollimatedBeam
from . import OpticalBench
from .lenses import *

    
class OpticalBench(object):
    """docstring for OpticalBench"""
    def __init__(self, x, y, z):
        super(OpticalBench, self).__init__()
        
        self.x = x
        self.y = y
        self.z = z
        
        self.boundary_list = self._makeBoundaries()
        self.source_list = []
        self.screen_list = []
        self.element_list = []
        
        self.render_limit = 5
    def _makeBoundaries(self):
        """docstring for makeBoundaries"""
        
        x = self.x
        y = self.y
        z = self.z
        
        return [Wall([2*x, y, z], [0,y,0],[0,0,z]),
                Wall([0, y, z], [0,y,0],[0,0,z]),
                Wall([x, 2*y, z], [x,0,0],[0,0,z]),
                Wall([x, 0, z], [x,0,0],[0,0,z]),
                Wall([x, y, 2*z], [x,0,0],[0,y,0]),
                Wall([x, y, 0], [x,0,0],[0,y,0]),]

    
    
    def Render(self):
        """Main render loop to perform ray tracing"""
        
        self.interactors = self.boundary_list + self.screen_list + self.element_list
        
        for source in self.source_list:
            for ray in source:
                i = 0
                while ray.isTerminated is False and i < self.render_limit:
                    ray = self._trace(ray)
                    i+=1
                    
    def _trace(self, ray):
        """Computes the distance along the ray to each interactor
            then propagates the ray through the one it hits first"""
        
        distances = []
        for obj in self.interactors:
            distances.append(obj.distance(ray))            
        distances  = np.array(distances)
        distances[np.logical_or(np.isnan(distances), distances < 0)] = float('inf')
        return self.interactors[np.argmin(distances)].propagate_ray(ray)
    
    
    def addSource(self, source):
        """Adds a source object to the bench"""
        assert isinstance(source, Source), "Argument must be a Source!"
        self.source_list.append(source)
    
    def addElement(self, element):
        """Adds an optical element object to the bench"""
        assert isinstance(element, OpticalElement), "Argument must be a Source!"
        self.element_list.append(element)
        
    def drawBench(self):
        """docstring for drawBench"""
        
        fig = plt.figure()
        ax = fig.add_axes([0.05,0.05, 0.9, 0.9])
        ax.set_xlim(-0.4*self.x,2.4*self.x)
        ax.set_ylim(-0.4*self.y,2.4*self.y)
        
        for b in self.boundary_list[:4]:
            b.drawBench(ax)
        
        for e in self.element_list[:4]:
            e.drawBench(ax)
        
        for s in self.source_list:
            for r in s:
                r.drawBench(ax)
glass = lambda l : 1.5046 + 4200./l**2
water =  lambda l : 1.319 + 6878./l**2        
air = lambda l : 1.                                
                
        
def test():
    ob = OpticalBench(0.5,0.5,.5)
    ob.addSource(CollimatedBeam([0.0, 0.5, 0.5], [1,0,0], 0.001, 35, 'white'))
    ob.addElement(BiConvex([0.2,0.5,0.5], 0.25, [-1.,0,0], 0.001, glass))
    ob.addElement(Screen([.85, 0.5, 0.5], [0,0.75,0],[0,0,0.5]))

    
    #for i in np.linspace(0,.99,30):
    #    ob.addSource(SingleRay([0.1, i, 0.5], [1,0,0], 500))
    #ob.addElement(PlaneRefraction([0.55, 0.5, 0.5], [0,-.2,0],[0,0,0.2], glass,air))
    #ob.addElement(Cube([0.55, 0.3, 0.5], [0.1, 0,0],[0,0.1,0],[0,0,0.1], water))
    
    #ob.addSource(CollimatedBeam([0.1, 0.1, 0.5], [1,.5,0], 0.001, 25, "white"))    
    #ob.addElement(CylindricalRefraction([0.5, 0.5, 0.5], L=[0.2,0,0], r=0.1, n1=air, n2=glass))    
    #ob.addElement(SphericalRefraction([0.5, 0.5, 0.5], [0.2,0,0],theta=np.pi, n1=air, n2=water))
    #ob.addElement(Mirror([0.5, 0.5, 0.5], [0,0.1,0],[0,0,0.2]))
    ##ob.addElement(Mirror([0.2, 0.7, 0.5], [0,0,.1],[0,.1,0]))
    
    
    #ob.addElement(PlaneRefraction([0.55, 0.6, 0.5], [0,0,0.2], [.3,0,0], air, glass))
    #ob.addElement(PlaneRefraction([0.55, 0.4, 0.5], [.3,0,0],[0,0,0.2], air,glass))
    #ob.addElement(PlanoConvex([0.5,0.5,0.5], 0.2, [-0.3,0,0], 0.05, glass))
    #ob.addElement(SphericalRefraction([1, 0.5, 0.5], [-0.3,0,0],theta=0.5, n1=air, n2=glass))
    
    #return ob
    ob.Render()
    ob.drawBench()
    return ob
    
def lensmakers(r1, r2, n, d):
    f = 1./( (n-1)*( 1./r1 - 1./r2 + (n-1)*d/(n*r1*r2)) )
    return f