from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from .utils import isVector, normalise, unifDisk, wavelengthToHex

air = lambda l: 1.


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
    def _reset(self):
        """docstring for _reset"""
        self.ray_list = self._initRays()
        
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
            self.u = normalise(np.array([1, 0, -self.direction[0] / self.direction[2]]))
        elif self.direction[1] != 0:
            self.u = normalise(np.array([1, -self.direction[0] / self.direction[1], 0]))
        elif self.direction[0] != 0:
            self.u = normalise(np.array([-self.direction[1] / self.direction[0], 1, 0]))
        else:
            print "invalid direction"
        self.v = normalise(np.cross(self.direction, self.u))

        self.ray_list = self._initRays()

    def _initRays(self):
        """docstring for _initRays"""
        pos = [self.centre + self.u * r * np.cos(t) + self.v * r * np.sin(t) for r, t in unifDisk(self.N, self.radius, 6)]
        if self.wavelength == 'white':
            return [Ray(p, self.direction, np.random.randint(400, 700)) for p in pos]
        else:
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

class PointSource(Source):

    """docstring for PointSource"""

    def __init__(self, centre, N,  wavelength):
        super(PointSource, self).__init__(centre)
        self.wavelength = wavelength
        self.N = N
        self.ray_list = self._initRays()

    def _initRays(self):
        """docstring for _initRays"""
        K = [[np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)] 
                for t in np.linspace(0,np.pi, self.N) 
                for p in np.linspace(0,2*np.pi, self.N)]
        return [Ray(self.centre, k, self.wavelength) for k in K] 


class ConicalSource(Source):

    """docstring for ConicalSource"""

    def __init__(self, centre, r,n, N,  wavelength):
        super(ConicalSource, self).__init__(centre)
        self.wavelength = wavelength
        self.r = r
        self.n = isVector(n)
        self.N = N
        self.ray_list = self._initRays()

    def _initRays(self):
        """docstring for _initRays"""
        K =[self.n + np.array([0, r * np.cos(t), r * np.sin(t)]) for r, t in unifDisk(self.N, self.r, 6)]
        if self.wavelength == 'white':
            return [Ray(self.centre, k, np.random.randint(400, 700)) for k in K]
        else:                        
            return [Ray(self.centre, k, self.wavelength) for k in K]
        
        
class Ray(object):

    """docstring for Ray"""

    def __init__(self, p, k, wavelength):
        super(Ray, self).__init__()
        self.k_history = []
        self.vertices = isVector(p).reshape((1, 3))
        self.k = isVector(k)        
        self.wavelength = wavelength
        self.isTerminated = False
        self._n = air(self.wavelength)  # Index of refraction of the medium currently in
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
        self.k_history.append(normalise(new))
        self._k = normalise(new)

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, new):
        self._n = new
        # print "ray.n = %s" %(self._n)

    def draw(self, ax=None):
        """docstring for draw"""
        if ax == None:
            fig = plt.figure()
            ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])

        ax.plot(self.vertices[:, 0], self.vertices[:, 1], '.-', color=wavelengthToHex(self.wavelength))

    def append(self, p, k):
        """Update the position and direction of the ray"""
        self.vertices = np.vstack((self.vertices, isVector(p)))
        self.k = isVector(k)
        # print "k = %s, p = %s" % (str(self.k), str(self.p))

    def __str__(self):
        s = "Ray:\nVertices: " + str(self.vertices) + "\nVector: " + str(self.k) + "\nMedium: " + str(self.n)
        return s
