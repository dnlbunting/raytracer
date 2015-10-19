#!/usr/bin/env python
from abc import ABCMeta, abstractmethod
from . import n_air
import numpy as np
import matplotlib.pyplot as plt
from .optical_elements import Wall, OpticalElement, Mirror, Screen, SphericalWall,Cube, SphericalMirror,PlaneInterface
from .sources import Source, SingleRay, CollimatedBeam
from .utils import rotPlane
    
    
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
        
        self.render_limit = 4
        
    def _makeBoundaries(self):
        """docstring for makeBoundaries"""
        
        x = self.x
        y = self.y
        z = self.z
        
        return [Wall([2*x, y, z], [0,y,0],[0,0,z]),
                Wall([0, y, z], [0,-y,0],[0,0,z]),
                Wall([x, 2*y, z], [x,0,0],[0,0,z]),
                Wall([x, 0, z], [-x,0,0],[0,0,z]),
                Wall([x, y, 2*z], [x,0,0],[0,y,0]),
                Wall([x, y, 0], [-x,0,0],[0,y,0]),]

    
    
    def Render(self):
        """Main render loop to perform ray tracing"""
        
        self.interactors = self.boundary_list + self.screen_list + self.element_list
        
        for source in self.source_list:
            for ray in source:
                i = 0
                self.last_interaction_idx = None
                while ray.isTerminated is False and i < self.render_limit:
                    ray = self._trace(ray)
                    i+=1
                           
    def _trace(self, ray):
        """Computes the distance along the ray to each interactor
            then propagates the ray through the one it hits first"""

        distances = []
        # Calculate distances to all objects before current
        for i in range(len(self.interactors)): 
            if i == self.last_interaction_idx:   
                # Distance to the element we interacted with last 
                # can never be 0 to avoid endless loops             
                d = self.interactors[i].distance(ray)
                d[d == 0] = float('inf')
                distances.append(np.min(d))
            else:
                distances.append(np.min(self.interactors[i].distance(ray)))
        distances = np.array(distances)        
            
        if ray.n != n_air:
            # The ray is currently inside some medium
            # handle getting the ray out in a single render loop
                        
            if np.sum(distances == 0) == 1:
                # The ray is currently on an interface between two objects
                print "Refracting into object"

                # The object we're about to move into
                new_obj_idx = int(np.where(distances == 0)[0])
                if self.interactors[self.last_interaction_idx].isTIR(ray, self.interactors[new_obj_idx].ref_index):
                    # Ray is TIR
                    ray = self.interactors[self.last_interaction_idx].TIR(ray)
                else:
                    # Ray propagates into new object
                    ray = self.interactors[new_obj_idx].propagate_ray(ray)                
                    self.last_interaction_idx = new_obj_idx
            
            elif np.sum(distances == 0) == 0:
                print "Refracting into air"
                # The ray is about to be emmited back into the air
                if self.interactors[self.last_interaction_idx].isTIR(ray, lambda x:n_air):
                    ray = self.interactors[self.last_interaction_idx].TIR(ray)                    
                else:
                    ray = self.interactors[self.last_interaction_idx].exitToAir(ray)
            
            else:
                # This is an error
                raise Exception("FATAL ERROR - 3 way interface encountered")
                        
            
        else:
            print "Interacts normally"
            ray = self.interactors[np.argmin(distances)].propagate_ray(ray)
            self.last_interaction_idx = np.argmin(distances)
        
        return ray
    
    
    def addSource(self, source):
        """Adds a source object to the bench"""
        assert isinstance(source, Source), "Argument must be a Source!"
        self.source_list.append(source)
    
    def addElement(self, element):
        """Adds an optical element object to the bench"""
        assert isinstance(element, OpticalElement), "Argument must be an OpticalElement!"
        self.element_list.append(element)
        
    def drawBench(self):
        """docstring for drawBench"""
        
        fig = plt.figure()
        ax = fig.add_axes([0.05,0.05, 0.9, 0.9])
        ax.set_xlim(-0.4*self.x,2.4*self.x)
        ax.set_ylim(-0.4*self.y,2.4*self.y)
        ax.set_aspect('equal')
        
        for b in self.boundary_list[:4]:
            b.drawBench(ax)
        
        for e in self.element_list[:4]:
            e.drawBench(ax)
        
        for s in self.source_list:
            for r in s:
                r.drawBench(ax)


glass = lambda l : 1.5046 + 4200/l**2
water =  lambda l : 1.319 + 6878/l**2        
                
def test():
    
    ob = OpticalBench(.5,.5,.5)
    for i in np.linspace(-.99,.99,10):
        ob.addSource(SingleRay([0.1, 0.5, 0.5], [1,i,0], 500))
    ob.addSource(SingleRay([0.15, 0.3, 0.5], [1,0.5,0], 500))
    #ob.addSource(SingleRay([0.9, 0.75, 0.5], [-1,-.5,0], 100))
    ob.addElement(Cube([0.3, 0.5, 0.5], [0.1, 0,0],[0,0.1,0],[0,0,0.1], glass))
    ob.addElement(Cube([0.5, 0.5, 0.5], [0.1, 0,0],[0,0.1,0],[0,0,0.1], lambda x:1))
    ob.addElement(Cube([0.3, 0.7, 0.5], [0.1, 0,0],[0,0.1,0],[0,0,0.1], lambda x:1.45))
    
    #ob.addSource(CollimatedBeam([0.1, 0.5, 0.5], [1,0,0], 0.05, 10, 100))    
    
    #ob.addElement(SphericalMirror([0.3, 0.5, 0.5], [0.5, 0,0],0.5))

    #ob.addElement(PlaneInterface([0.5, 0.5, 0.5], [0,-0.1,0],[0,0,0.2], 1.3))
    #ob.addElement(Mirror([0.5, 0.5, 0.5], *rotPlane(0.1,0.1,np.pi/4.)))
    
    #ob.addElement(Screen([0.5, 0.9, 0.5], [0.5,0,0],[0,0,0.5]))
    
    ob.Render()
    ob.drawBench()
    return ob
    

