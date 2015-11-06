#!/usr/bin/env python
from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt

from .optical_elements import *
from .sources import Source, SingleRay, CollimatedBeam
from .lenses import *
from .utils import validDistanceArray

eps = np.finfo(np.float64).eps


class OpticalBench(object):

    """OpticalBench is the main class of the raytracer module. To create a simulation generally
    create an OpticalBench object and add sources/optical elements to it then call Render() and draw()"""

    def __init__(self, x, y, z, verbose=False, render_limit=20):
        '''
        Creates an OpticalBench object to which optical elements can be added to build up 
        a simulation scene
        
        x,y,z -  The dimensions of the scene
        verbose - Print each ray during the render loop - default:False
        render_limit - Maximum number of times to trace a single ray. Used to prevent 
                       infinite loops ie two mirrors facing each other. Throws a RenderLimitExceeded
                       exception if violated
                    
        
        '''
        super(OpticalBench, self).__init__()

        self.x = x
        self.y = y
        self.z = z

        self.verbose = verbose

        self.boundary_list = self._makeBoundaries()
        self.source_list = []
        self.screen_list = []
        self.element_list = []

        self.render_limit = render_limit

    def _makeBoundaries(self):
        """Private function. Builds the edges of the simulation area"""

        x = self.x
        y = self.y
        z = self.z

        return [Wall([2 * x, y, z], [0, y, 0], [0, 0, z]),
                Wall([0, y, z], [0, y, 0], [0, 0, z]),
                Wall([x, 2 * y, z], [x, 0, 0], [0, 0, z]),
                Wall([x, 0, z], [x, 0, 0], [0, 0, z]),
                Wall([x, y, 2 * z], [x, 0, 0], [0, y, 0]),
                Wall([x, y, 0], [x, 0, 0], [0, y, 0]), ]

    def Render(self):
        """Main render loop to perform ray tracing"""

        self.interactors = self.boundary_list + self.screen_list + self.element_list

        for source in self.source_list:
            source._reset()
            for ray in source:
                i = 0
                while ray.isTerminated is False:
                    if i > self.render_limit:
                        raise Exception("Render limit exceeded")
                    if self.verbose:
                        print "\n"
                        print "Loop " + str(i)
                        print ray
                    ray = self._trace(ray)
                    i += 1
                if self.verbose:
                    print "*" * 30

    def _trace(self, ray):
        """Private. Computes the distance along the ray to each interactor
            then propagates the ray through the one it hits first"""

        distances = []
        for obj in self.interactors:
            distances.append(obj.distance(ray))
        distances = np.array(distances)
        if self.verbose:
            print "Distances: " + str(distances)
        # Assert that distances must be >= 0 and not nan
        distances = validDistanceArray(distances)
        return self.interactors[np.argmin(distances)].propagate_ray(ray)

    def add(self, element):
        """Adds an optical element object to the bench"""
        if isinstance(element, Source):
            self.source_list.append(element)
        elif isinstance(element, Screen):
            self.screen_list.append(element)
        elif isinstance(element, OpticalElement):
            self.element_list.append(element)
        else:
            raise Exception("Object must be subclassed from OpticalElement or Source")

    def draw(self, ax=None):
        """Draws a xy projection of the OpticalBench. Generally call Render() first to
            populate with rays, although can be called before render to check the placement of objects etc"""
        if ax == None:
            fig = plt.figure()
            ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
            ax.set_xlim(-0.4 * self.x, 2.4 * self.x)
            ax.set_ylim(-0.4 * self.y, 2.4 * self.y)
            ax.set_aspect('equal')

        for b in self.boundary_list[:4]:
            b.draw(ax)

        for e in self.element_list[:4]:
            e.draw(ax)

        for e in self.screen_list:
            e.draw(ax)
            e.fig.show()

        for s in self.source_list:
            for r in s:
                r.draw(ax)

        return ax


