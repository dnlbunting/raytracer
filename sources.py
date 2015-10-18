from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from .utils import isVector, normalise, unifDisk
from . import n_air



class Source(object):
    """docstring for Source"""
    __metaclass__ = ABCMeta
    
    def __init__(self, centre):
        super(Source, self).__init__()
        
        self.centre = centre
        self.ray_list = []
        
    def __getitem__(self, key):
        """Allows for direct access and iteration over the rays emitted from this source"""
        return self.ray_list[key]
    
    @abstractmethod
    def _initRays(self):
        pass    

class CollimatedBeam(Source):
    """docstring for CollimatedBeam"""
    def __init__(self, centre, direction, radius, N, wavelength):
        super(CollimatedBeam, self).__init__(centre)
        
        self.centre = centre
        self.direction = normalise(direction)
        self.radius = radius
        self.N = N
        self.wavelength = wavelength
        
        
        # Create an orthonormal basis on the disk
        if self.direction[2] != 0:
            self.u = normalise(np.array([1,0, -self.direction[0]/self.direction[2]]))
        elif self.direction[1] != 0:
            self.u = normalise(np.array([1, -self.direction[0]/self.direction[1],0]))
        elif self.direction[0] != 0:
            self.u = normalise(np.array([-self.direction[1]/self.direction[0],1,0]))
        else:
            print "invalid direction" 
        self.v = normalise(np.cross(self.direction, self.u))
        
        self.ray_list =  self._initRays()
    def _initRays(self):
        """docstring for _initRays"""
        pos = [self.centre + self.u * r*np.cos(t) + self.v *r*np.sin(t) for r,t in unifDisk(10,self.radius,6)]
        return [Ray(p, self.direction, self.wavelength) for p in pos] 
        
    
class SingleRay(Source):
    """docstring for SingleRay"""
    def __init__(self, centre, direction,  wavelength):
        super(SingleRay, self).__init__(centre)
        self.direction = isVector(direction)
        self.wavelength = wavelength
        
        self.ray_list = self._initRays()
    def _initRays(self):
        """docstring for _initRays"""
        return([Ray(self.centre, self.direction, self.wavelength)])

class Ray(object):
    """docstring for Ray"""
    
    def __init__(self, p, k,wavelength):
        super(Ray, self).__init__()
        
        self.vertices = isVector(p).reshape((1,3))
        self.k = isVector(k)
        self.wavelength = wavelength
        self.isTerminated = False
        self.n = n_air # Index of refraction of the medium currently in
                   # TODO make less of hack
        
    
    @property
    def p(self):
        """The current point of the ray, i.e. most recent vertex"""
        return self.vertices[-1]
    
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, new):
        self._k = normalise(new)
                        
    def drawBench(self, ax=None):
        """docstring for drawBench"""
        if ax == None:
            fig = plt.figure()
            ax = fig.add_axes([0.05,0.05, 0.9, 0.9])
            
        ax.plot(self.vertices[:,0], self.vertices[:,1],'.-')
    
    def append(self, p, k):
        """Update the position and direction of the ray"""
        self.vertices = np.vstack((self.vertices, isVector(p)))
        self.k = isVector(k)