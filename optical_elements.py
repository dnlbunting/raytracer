from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D

from .utils import isVector, normalise, wavelengthToHex, validDistanceArray
from .interactions import RefractionMixin, ReflectionMixin, AbsorptionMixin

air = lambda l: 1.
eps = np.finfo(np.float64).eps

# TODO Add kwargs to draw that pass to Polygon to allow for more complex rendering ie colors for materials

# Abstact bases


class OpticalElement(object):

    """Abstract base class defining the common interface that all optical elements
        must implement and subclass.

        Also defines draw() which should be called first in all derived draw
        implementations to handle setting up new axes if none are provided"""
    __metaclass__ = ABCMeta

    def __init__(self, centre):
        super(OpticalElement, self).__init__()

        self.centre = isVector(centre)

    @abstractmethod
    def distance(self):
        pass

    @abstractmethod
    def propagate_ray(self):
        pass

    @abstractmethod
    def draw(self, ax=None):
        if ax == None:
            fig = plt.figure()
            ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], aspect='auto')
            ax.set_xlim(-2, 2)
            ax.set_ylim(-2, 2)

        return ax


class CompositeObject(OpticalElement):

    """docstring for CompositeObject"""

    def __init__(self, centre):
        super(CompositeObject, self).__init__(centre)
        self.components = []

    def add(self, comp):
        """Add a new component to the composite object"""
        if not isinstance(comp, OpticalElement):
            raise Exception("Component must be an optical element!")

        self.components.append(comp)
        
    def distance(self, ray):
        """Distance of the composite is the min distance of its components"""
        d = np.array([s.distance(ray) for s in self.components])
        d = validDistanceArray(d)
        return np.min(d)

    def propagate_ray(self, ray):
        """Pass on propagate_ray to the component closest to the ray"""
        d = np.array([s.distance(ray) for s in self.components])
        d = validDistanceArray(d)
        return self.components[np.argmin(d)].propagate_ray(ray)

    def draw(self, ax=None):
        """docstring for draw"""
        ax = super(CompositeObject, self).draw(ax)
        for s in self.components:
            s.draw(ax)
        return ax


# Cyclindrical Objects
class Cylidrical(OpticalElement):

    """docstring for Cylidrical"""

    def __init__(self, centre, L, r):
        super(Cylidrical, self).__init__(centre)
        self.L = isVector(L)
        if L[2] != 0:
            raise Exception("Only cylinders in the xy plane are currently supported :(")
        self.r = r
        self.R = normalise(np.cross([0, 0, 1], self.L))
        self.points = [self.centre + self.L + self.r * self.R,
                       self.centre + self.L - self.r * self.R,
                       self.centre - self.L - self.r * self.R,
                       self.centre - self.L + self.r * self.R]

    def draw(self, ax=None):
        ax = super(Cylidrical, self).draw(ax)
        poly = Polygon([x[:2] for x in self.points], True, facecolor='none')
        ax.add_patch(poly)

    def distance(self, ray):

        r_prime = self.centre - ray.p

        # Projection onto cylinder cross section
        if np.abs(np.dot(ray.k, normalise(self.L))) == 1.:
            # Ray is parrallel to the cylinder
            return float('inf')

        x = np.dot(r_prime, [0, 0, 1]) * np.array([0, 0, 1]) + np.dot(r_prime, self.R) * self.R
        k = normalise(np.dot(ray.k, [0, 0, 1]) * np.array([0, 0, 1]) + np.dot(ray.k, self.R) * self.R)

        # If the sqrt is imaginary it never intersects
        sq = np.dot(x, k)**2 - (np.linalg.norm(x)**2 - self.r**2)
        if sq > 0:
            l1 = np.dot(x, k) + np.sqrt(sq)
            l2 = np.dot(x, k) - np.sqrt(sq)
        else:
            return float('inf')

        # Try the closest intersection first
        d = sorted(np.array([l1, l2]) / np.dot(ray.k, k))
        d = validDistanceArray(d)
        if d[0] != float('inf') and np.abs(np.dot(self.centre - ray.p + d[0] * ray.k, self.L) / np.linalg.norm(self.L)**2) <= 1:
            return d[0]

        # Now try the further away one
        elif d[1] != float('inf') and np.abs(np.dot(self.centre - ray.p + d[1] * ray.k, self.L) / np.linalg.norm(self.L)**2) <= 1:
            return d[1]
        else:
            return float('inf')


class CylindricalRefraction(Cylidrical, RefractionMixin):

    """docstring for CylindricalRefraction"""

    def __init__(self, centre, L, r, n1, n2):
        super(CylindricalRefraction, self).__init__(centre, L, r)
        self.n1 = n1
        self.n2 = n2

    def propagate_ray(self, ray):
        print "WARNING- Cyclindrical refraction"
        d = self.distance(ray)
        p = ray.p + d * ray.k
        z = np.cross(p - self.centre, self.L)
        n = normalise(np.cross(self.L, z))

        return self._refract(ray, n)


class CylindricalWall(Cylidrical, AbsorptionMixin):

    """docstring for CylindricalWall"""

    def __init__(self, centre, L, r):
        super(CylindricalWall, self).__init__(centre, L, r)

    def propagate_ray(self, ray):
        return self._absorb(ray)


# Spherical Objects
class Spherical(OpticalElement):

    """docstring for Spherical"""

    def __init__(self, centre, R, theta=None, height=None):
        super(Spherical, self).__init__(centre)

        # scalar radius
        self.r = np.linalg.norm(R)
        # Unit vector defining radial direction
        self.R = normalise(R)

        self.theta = theta
        self.height = height
        self.t_0 = np.arctan2(self.R[1], self.R[0])

        if (theta is None and height is None) or (theta is not None and height is not None):
            raise Exception("Supply exactly one of theta or height")
        elif theta is None:
            self.theta = np.arcsin(self.height / self.r)
        elif height is None:
            if self.theta < 0.5 * np.pi:
                self.height = self.r * np.sin(self.theta)
            else:
                self.height = self.r
        assert np.logical_and(self.theta > 0, theta <= np.pi), "0 < theta <= pi"

    def draw(self, ax=None):
        """docstring for draw"""
        ax = super(Spherical, self).draw(ax)

        if self.R[2] == 0:
            # If the spherical isn't tilted in the z direction we can
            # use a simplified plotting algorithm since we know the boundary

            T = np.linspace(self.t_0 - self.theta, self.t_0 + self.theta, 100)
            points = [[self.r * np.cos(t) + self.centre[0], self.r * np.sin(t) + self.centre[1]] for t in T]
            poly = Polygon(points, True, facecolor='none')
            ax.add_patch(poly)

        # TODO Make it so that the spherical can have a z comp

        else:
            raise NotImplementedError()

    def plot3D(self):
        """docstring for plot3D"""
        raise NotImplementedError()

    def distance(self, ray):
        r = ray.p - self.centre

        # If the sqrt is imaginary it never intersects
        sq = np.dot(r, ray.k)**2 - (np.linalg.norm(r)**2 - self.r**2)
        if sq >= 0:
            l1 = -np.dot(r, ray.k) + np.sqrt(sq)
            l2 = -np.dot(r, ray.k) - np.sqrt(sq)
        else:
            return float('inf')

        # Try the closest intersection first
        d = sorted([l1, l2])
        d = validDistanceArray(d)
        if d[0] != float('inf') and self.isIntersect(ray.p + d[0] * ray.k - self.centre):
            return d[0]
        # Now try the further away one
        elif d[1] != float('inf') and self.isIntersect(ray.p + d[1] * ray.k - self.centre):
            return d[1]
        else:
            return float('inf')

    def isIntersect(self, p):
        """docstring for isIntersect"""

        ang = np.arccos(np.dot(p, self.R) / np.linalg.norm(p))
        z_cond = np.abs(p[2]) <= self.height
        return np.logical_and(z_cond, np.abs(ang) <= self.theta)


class SphericalWall(Spherical, AbsorptionMixin):

    """A spherical with total absorption"""

    def __init__(self, centre,  R, theta=None, height=None):
        super(SphericalWall, self).__init__(centre, R, theta, height)

    def propagate_ray(self, ray):
        return self._absorb(ray)


class SphericalMirror(Spherical, ReflectionMixin):

    """A totally reflective spherical (both sides)"""

    def __init__(self, centre,  R, theta=None, height=None):
        super(SphericalMirror, self).__init__(centre, R, theta, height)

    def propagate_ray(self, ray):
        """Implements propagate_ray for Wall by terminating ray"""

        d = self.distance(ray)
        p = ray.p + d * ray.k

        # Normal at the point of intersection
        n = normalise(p - self.centre)

        return self._reflect(ray, n)


class SphericalRefraction(Spherical, RefractionMixin):

    """A spherical with medium on one side and with air on the other"""

    def __init__(self, centre,  R, n1, n2, theta=None, height=None):
        super(SphericalRefraction, self).__init__(centre, R, theta, height)
        self.n1 = n1
        self.n2 = n2

    def propagate_ray(self, ray):
        """Implements propagate_ray for SphericalRefraction by calculating the
            refracted wave vector using the formula from wikipedia"""

        d = self.distance(ray)
        n = normalise(ray.p + d * ray.k - self.centre)
        # print n
        return self._refract(ray, n)


# Planar Objects
class Plane(OpticalElement):

    """Plane optical element abstract class, used for OpticalBench
        boundaries, mirrors etc

        pos - Vector to the center of the plane
        a,b - Vectors describing the distance from the center to each edge,
              *must be orthogonal*. Side lengths are then 2a, 2b"""

    def __init__(self, centre,  a, b):
        super(Plane, self).__init__(centre)

        assert np.dot(a, b) < eps, "Vectors must be orthogonal"

        self.a = isVector(a)
        self.b = isVector(b)

        self.normal = normalise(np.cross(a, b))

        # The four corners of the plane
        self.points = [self.centre + a + b,
                       self.centre - a + b,
                       self.centre - a - b,
                       self.centre + a - b]

    def draw(self, ax=None):
        """Draws a projection of the object on the xy plane"""
        ax = super(Plane, self).draw(ax)
        poly = Polygon([x[:2] for x in self.points], True, facecolor='none')
        ax.add_patch(poly)
    
    def distance(self, ray):
        """Distance from ray to self along the ray's path"""

        r = self.centre - ray.p
        top = np.dot(self.normal, r)
        bottom = np.dot(ray.k, self.normal)
        if np.abs(top) < eps or np.abs(bottom) < eps:
            return float('inf')

        d = top / bottom
        p = ray.p + d * ray.k
        d = d if self.isIntersect(p) else float('inf')
        return d

    def isIntersect(self, p):
        """Decides if p is in the plane or not by projecting the intersection point
            on the infinite plane onto the a,b axis than span the finite plane
            and seeing if either component is >1"""
        alpha = np.dot(self.a, (p - self.centre)) / np.linalg.norm(self.a)**2
        beta = np.dot(self.b, (p - self.centre)) / np.linalg.norm(self.b)**2

        return np.logical_and(np.abs(alpha) <= 1, np.abs(beta) <= 1)


class PlaneRefraction(Plane, RefractionMixin):

    """A plane with some medium with refractive index n1 on one side
         and n2 on the other. The normal vector is defined to point towards
         the medium with n1"""

    def __init__(self, centre,  a, b, n1, n2):
        super(PlaneRefraction, self).__init__(centre,  a, b)
        self.n1 = n1
        self.n2 = n2

    def propagate_ray(self, ray):
        """Implements propagate_ray for PlaneRefraction by calculating the
            refracted wave vector using the formula from wikipedia"""
        return self._refract(ray, self.normal)


class Mirror(Plane, ReflectionMixin):

    """A totally reflective plane (both sides)"""

    def __init__(self, centre,  a, b):
        super(Mirror, self).__init__(centre,  a, b)

    def propagate_ray(self, ray):
        """Implements propagate_ray for Mirror by reflecting the wavevector of the ray
            about the normal of the plane.
            """

        return self._reflect(ray, self.normal)


class Wall(Plane, AbsorptionMixin):

    """A plane with total absorption"""

    def __init__(self, centre,  a, b):
        super(Wall, self).__init__(centre,  a, b)

    def propagate_ray(self, ray):
        """Implements propagate_ray for Wall by terminating ray"""
        return self._absorb(ray)


class Cube(CompositeObject):

    """docstring for Cube"""

    def __init__(self, centre,  a, b, c, ref_index):
        super(Cube, self).__init__(centre)
        self.a = isVector(a)
        self.b = isVector(b)
        self.c = isVector(c)

        self.ref_index = ref_index

        self.side_centres = [self.centre + self.c,
                             self.centre - self.c,
                             self.centre + self.b,
                             self.centre - self.b,
                             self.centre + self.a,
                             self.centre - self.a]

        sides = [PlaneRefraction(self.side_centres[0],  self.a, self.b, air, self.ref_index),
                 PlaneRefraction(self.side_centres[1],  -self.a, self.b, air, self.ref_index),
                 PlaneRefraction(self.side_centres[2],  -self.a, self.c, air, self.ref_index),
                 PlaneRefraction(self.side_centres[3],  self.a, self.c, air, self.ref_index),
                 PlaneRefraction(self.side_centres[4],  self.b, self.c, air, self.ref_index),
                 PlaneRefraction(self.side_centres[5],  -self.b, self.c, air, self.ref_index)]
        for s in sides:
            self.add(s)


class Screen(Plane):

    """An output screen, records all incident beams as pixels and renders them"""

    def __init__(self, centre,  a, b, term = False):
        super(Screen, self).__init__(centre,  a, b)
        self.term = term
        self.pixels = []

    def propagate_ray(self, ray):
        """Implements propagate_ray for Screen, projects the ray's intersection with
            the screen onto plane orthogonal coordinates and draws it as a point on the
            associated mpl figure"""
        d = self.distance(ray)
        ray.append(ray.p + d * ray.k, ray.k)
        if self.term:
            ray.isTerminated = True

        alpha = np.dot(self.a, (ray.p - self.centre)) / np.linalg.norm(self.a)
        beta = np.dot(self.b, (ray.p - self.centre)) / np.linalg.norm(self.b)
        c = wavelengthToHex(ray.wavelength)

        self.pixels.append([alpha, beta, c])

        return ray

    def draw(self, ax=None):
        """Draws a projection of the object on the xy plane"""
        ax = super(Screen, self).draw(ax)
        
        self.fig = plt.figure()
        self.ax = self.fig.add_axes([0.05, 0.05, 0.9, 0.9], aspect='auto')
        A = np.linalg.norm(self.a)
        B = np.linalg.norm(self.b)
        self.ax.set_xlim(-1 * A, A)
        self.ax.set_ylim(-1 * B, B)
        self.ax.set_aspect('equal')
        
        for p in self.pixels:
            self.ax.scatter(p[0], p[1], c=p[2], marker='.', lw=0)


class Mask(Plane):

    """"""

    def __init__(self, centre,  a, b, mask):
        super(Mask, self).__init__(centre,  a, b)

        self.pixels = []
        self.mask = np.load(mask)
        self.w = 0.5*self.mask.shape[0]
        self.h = 0.5*self.mask.shape[1]
        self.border = 100
    
    def propagate_ray(self, ray):
        """"""
        d = self.distance(ray)
        ray.append(ray.p + d * ray.k, ray.k)

        # Plane orthogonal coordinates, rescale into pixels
        alpha = np.dot(self.a, (ray.p - self.centre)) / np.linalg.norm(self.a)**2
        beta = np.dot(self.b, (ray.p - self.centre)) / np.linalg.norm(self.b)**2
        
        alpha = (self.border+self.w) * 2*alpha
        beta = (self.border+self.h) * 2*beta
        if np.abs(alpha) > self.w or np.abs(beta) > self.h:
            # Outside of image area, terminate
            ray.isTerminated = True
            return ray

        elif self.mask[self.w + np.floor(alpha), self.h+np.floor(beta)] > 0 :
            ray.isTerminated = False
            return ray
        else :
            ray.isTerminated = True
            return ray
        

        return ray


