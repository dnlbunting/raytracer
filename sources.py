from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from .utils import isVector, normalise, unifDisk, wavelengthToHex

air = lambda l: 1.


class Source(object):

    """ABC defining common interface for all source like objects
        
        ray_list - stores all rays associated with this source, exposed for direct access/iteration
        _initRays - populates ray_list 
        _reset - re-populates ray_list """
        
    __metaclass__ = ABCMeta

    def __init__(self, centre):
        super(Source, self).__init__()

        self.centre = centre
        self.ray_list = []

    def __getitem__(self, key):
        """Allows for direct access and iteration over the rays emitted from this source"""
        return self.ray_list[key]
    
    def _reset(self):
        """Returns source to its original state ie before call to Render()"""
        self.ray_list = self._initRays()
        
    @abstractmethod
    def _initRays(self):
        pass


class CollimatedBeam(Source):

    """A beam of collimated light
    
        centre -  centre point of the beam 
        direction -  direction of propagation vector
        radius  - radius of the beam 
        N - number of rings
        wavelength -  wavelength in nm, can be 'white' to produce beam with random colored rays """

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
        pos = [self.centre + self.u * r * np.cos(t) + self.v * r * np.sin(t) for r, t in unifDisk(self.N, self.radius, 6)]
        if self.wavelength == 'white':
            return [Ray(p, self.direction, np.random.randint(400, 700)) for p in pos]
        else:
            return [Ray(p, self.direction, self.wavelength) for p in pos]
    



class SingleRay(Source):

    """Wrapper for a  single ray object
    
       centre -  centre point of the beam 
       direction -  direction of propagation vector
       wavelength -  wavelength in nm"""

    def __init__(self, centre, direction,  wavelength):
        super(SingleRay, self).__init__(centre)
        self.direction = isVector(direction)
        self.wavelength = wavelength

        self.ray_list = self._initRays()

    def _initRays(self):
        """docstring for _initRays"""
        return([Ray(self.centre, self.direction, self.wavelength)])

class PointSource(Source):

    """Point source, emits rays uniformly in all directions
       
       centre -  centre point of the beam 
       direction -  direction of propagation vector
       N - angular number density of rays (rays/pi)
       wavelength -  wavelength in nm"""

    def __init__(self, centre, N,  wavelength):
        super(PointSource, self).__init__(centre)
        self.wavelength = wavelength
        self.N = N
        self.ray_list = self._initRays()

    def _initRays(self):
        K = [[np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)] 
                for t in np.linspace(0,np.pi, self.N) 
                for p in np.linspace(0,2*np.pi, 2*self.N)]
        return [Ray(self.centre, k, self.wavelength) for k in K] 


class ConicalSource(Source):

    """ConicalSource, emits rays from centre point so that they uniformly cover 
        a disk of radius r,  1 unit away from centre in the direction n
        
        centre -  centre point of the beam 
        r - radius of the disk the rays cover 1 unit away from the centre
        direction -  direction of the centre propagation vector
        N -  number  of rays
        wavelength -  wavelength in nm"""

    def __init__(self, centre, r, direction, N,  wavelength):
        super(ConicalSource, self).__init__(centre)
        self.wavelength = wavelength
        self.r = r
        self.direction = isVector(direction)
        self.N = N
        self.ray_list = self._initRays()

    def _initRays(self):
        """docstring for _initRays"""
        K =[self.direction + np.array([0, r * np.cos(t), r * np.sin(t)]) for r, t in unifDisk(self.N, self.r, 6)]
        if self.wavelength == 'white':
            return [Ray(self.centre, k, np.random.randint(400, 700)) for k in K]
        else:                        
            return [Ray(self.centre, k, self.wavelength) for k in K]
        
        
class Ray(object):

    """Ray object"""

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
        """The current direction vector of the ray"""
        return self._k

    @k.setter
    def k(self, new):
        self.k_history.append(normalise(new))
        self._k = normalise(new)

    @property
    def n(self):
        """The refractive index of the medium the ray is currently travelling in"""
        return self._n

    @n.setter
    def n(self, new):
        self._n = new
        # print "ray.n = %s" %(self._n)

    def draw(self, ax=None):
        """Draws the xy projection of the rays path"""
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
