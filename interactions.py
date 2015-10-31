#!/usr/bin/env python
from .utils import isVector, normalise, wavelengthToHex
import numpy as np
air = lambda l: 1.
eps = np.finfo(np.float64).eps

"""This module contains Mixin classes designed to be subclassed by objects alongside one of
   the geometry ABCs in optical_elements in order to provide optical interactions requiring 
    only the incident normal."""


class ReflectionMixin(object):
    """ReflectionMixin provides geometry indenpendent reflection"""
    def _reflect(self, ray, n):
        '''Reflects an incident ray, ray in the plane with normal n'''
        d = self.distance(ray)
        p = ray.p + d * ray.k
        k_prime = ray.k - 2 * np.dot(ray.k, n) * n
        ray.append(p, k_prime)
        ray.isTerminated = False
        return ray


class RefractionMixin(object):

    """RefractionMixin provide geometry independent refraction"""

    def __init__(self):
        super(RefractionMixin, self).__init__()

    def _refract(self, ray, n):
        """Refracts an incident ray with a plane with normal n.
            Requires that the object subclassing this mixin has members n1 and n2, where n is
            taken to point from n2 -> n1.
            Direction of the refraction is decided automatically based on the angle between ray.k and n, 
            and we assert that this is in agreement with the medium information carried by the ray object """
        
        d = self.distance(ray)
        c = -np.dot(n, ray.k)

        # Set up the oreintation of the interface
        # print "Normal = %s" %(str(n))
        if c > 0:
            # print "Refracting n1 -> n2"
            # Ray is propagating n1 -> n2
            if np.abs(ray.n - self.n1(ray.wavelength)) > eps:
                raise Exception("Ray current refractive index %f does not match that of the interface %f" % (ray.n, self.n1(ray.wavelength)))
            r = self.n1(ray.wavelength) / self.n2(ray.wavelength)
            # print "r = %g / %g " % (self.n1(ray.wavelength), self.n2(ray.wavelength) )

        elif c < 0:
            # print "Refracting n2 -> n1"
            # Ray is propagating n2 -> n1
            if np.abs(ray.n - self.n2(ray.wavelength)) > eps:
                raise Exception("Ray current refractive index %f does not match that of the interface %f" % (ray.n, self.n2(ray.wavelength)))

            r = self.n2(ray.wavelength) / self.n1(ray.wavelength)
            # print "r = %g / %g " % (self.n2(ray.wavelength) ,self.n1(ray.wavelength))
            c = -c
            n = -n

        else:
            raise Exception("A fatal error has occurred - Ray is parallel to the interface")

        if 1 - (r**2) * (1 - c**2) < 0:
            # Ray is TIR
            k_prime = ray.k - 2 * np.dot(ray.k, n) * n
            ray.append(ray.p + d * ray.k, k_prime)
        else:
            # Ray is refracted
            k_prime = r * ray.k + (r * c - np.sqrt(1 - (r**2) * (1 - c**2))) * n
            ray.append(ray.p + d * ray.k, k_prime)
            # print "Medium changed from %g to %g" %(ray.n,ray.n/r )
            ray.n = ray.n / r

        ray.isTerminated = False

        return ray


class AbsorptionMixin(object):
    """Terminates rays"""
    def _absorb(self, ray):
        d = self.distance(ray)
        ray.append(ray.p + d * ray.k, ray.k)
        ray.isTerminated = True
        return ray
