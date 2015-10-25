#!/usr/bin/env python
from .utils import isVector, normalise, wavelengthToHex
import numpy as np
air = lambda l: 1.
eps = np.finfo(np.float64).eps


class ReflectionMixin(object):

    def _reflect(self, ray, n):

        d = self.distance(ray)
        p = ray.p + d * ray.k
        k_prime = ray.k - 2 * np.dot(ray.k, n) * n
        ray.append(p, k_prime)
        ray.isTerminated = False
        return ray


class RefractionMixin(object):

    """docstring for RefractionMixin"""

    def __init__(self):
        super(RefractionMixin, self).__init__()

    def _refract(self, ray, n):
        """The actual refraction/TIR is performed here as once the normal n has been
            defined by the geometry specific code the process is generic"""
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

    def _absorb(self, ray):
        """docstring for _absorb"""
        d = self.distance(ray)
        ray.append(ray.p + d * ray.k, ray.k)
        ray.isTerminated = True
        return ray
