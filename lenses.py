#!/usr/bin/env python
from .utils import isVector, normalise, wavelengthToHex
import numpy as np
from .optical_elements import *

air = lambda l: 1.
'''
Lens objects.

Note convex lenses have been tested way more than the concave, which are still unstable
'''
def lensmakers(r1, r2, n, d):
    f = 1. / ((n - 1) * (1. / r1 - 1. / r2 + (n - 1) * d / (n * r1 * r2)))
    return f

class PlanoConvex(CompositeObject):

    """PlanoConvex lens object
        
        Params:
        centre - Centre of the object
        height -  height of the lens
        R - radius of curvature vector
        thickness - thickness of the lens, excluding the thickness due to the curved surface
        ref_index - refractive index as function of wavelength. For constant use eg lambda l : 1.5 
        
        Members:
        f(wl) - The focal length at wavelength wl
        bf1/bfd2 (wl) -  The back focal distances. Ie distance from the edge of the lens to the focal point.
                        using bfd1 or 2 depends on whether the curved or planar surface faces the source.
        """

    def __init__(self, centre, height, R, thickness, ref_index):
        super(PlanoConvex, self).__init__(centre)
        self.height = height
        self.radius = np.linalg.norm(R)
        self.R = normalise(isVector(R))
        self.thickness = thickness
        self.ref_index = ref_index
        self.w = self.radius - np.sqrt(self.radius**2 - self.height**2)
        self.d = self.thickness + self.w
        
        # Build the composite
        self.add(SphericalRefraction(self.centre - (self.radius - self.w - 0.5 * self.thickness) * self.R,
                                     self.radius * self.R, height=self.height, n1=air, n2=ref_index))
        if self.thickness > 0:
            self.add(CylindricalRefraction(self.centre, self.R * self.thickness * 0.5,
                                       self.height, self.ref_index, air))

        a = self.height * np.array([0, 0, 1])
        b = self.height * normalise(np.cross(a, self.R))
        self.add(PlaneRefraction(self.centre - self.R * 0.5 * self.thickness, a, b, n1=air, n2=ref_index))

        # Analytic calculation of focal point
        self.f = lambda wl: lensmakers(float('inf'), -1. * self.radius, self.ref_index(wl), self.thickness + self.w)
        self.bfd1 = lambda wl : self.f(wl) * (1 - (self.ref_index(wl)-1)*self.d/(self.ref_index(wl)*self.radius) )
        self.bfd2 = self.f

class PlanoConcave(CompositeObject):

    """PlanoConcave lens object- Untested
        
        Params:
        centre - Centre of the object
        height -  height of the lens
        R - radius of curvature vector
        thickness - thickness of the lens, excluding the thickness due to the curved surface
        ref_index - refractive index as function of wavelength. For constant use eg lambda l : 1.5"""
        
    def __init__(self, centre, height, R, thickness, ref_index):
        super(PlanoConcave, self).__init__(centre)
        self.height = height
        self.radius = np.linalg.norm(R)
        self.R = normalise(isVector(R))
        self.ref_index = ref_index
        self.w = self.radius - np.sqrt(self.radius**2 - self.height**2)
        self.thickness = thickness + self.w

        self.add(SphericalRefraction(self.centre + (self.radius - self.w + 0.5 * self.thickness) * self.R,
                                     -self.radius * self.R, height=self.height, n1=ref_index, n2=air, draw_kwargs={'closed':False, 'facecolor':'none', 'edgecolor':'red'}))
        if self.thickness > 0:
            self.add(CylindricalRefraction(self.centre, self.R * self.thickness * 0.5,
                                           self.height, air, self.ref_index))

        a = self.height * np.array([0, 0, 1])
        b = self.height * normalise(np.cross(a, self.R))
        self.add(PlaneRefraction(self.centre - self.R * 0.5 * self.thickness, a, b, n1=air, n2=ref_index))


        self.f = lambda wl : lensmakers(float('inf'), 1. * self.radius, self.ref_index(wl), self.thickness + self.w)


class BiConvex(CompositeObject):

    """BiConvex lens. Both r1 and r2 are assumed to be positive and the surface with r1 is facing the source
    
        Params:
        centre - Centre of the object
        height -  height of the lens
        r1 - Scalar radius of curvature of the first surface (>0)
        r2 - Scalar radius of curvature of the second surface (>0)        
        R - Vector defining the optical axis of the lens
        thickness - thickness of the lens, excluding the thickness due to the curved surfaces
        ref_index - refractive index as function of wavelength. For constant use eg lambda l : 1.5
        
        Members:
        f(wl) - The focal length at wavelength wl
        bfd(wl) -  The back focal distance. Ie distance from the edge of the lens to the focal point.
                    This assumes the light is incident of side r1 so bfd is the distance from the edge of r2
                    to the focal length. The distance from the centre is then bfd+w2+0.5*thickness
        w1,w2 - The widths of the curved surfaces
        d - The total width of the lens
        """

    def __init__(self, centre, height, r1, r2, R, thickness, ref_index):
        super(BiConvex, self).__init__(centre)
        self.height = height
        
        self.radius1 = np.abs(r1)
        self.R1 = -1.*normalise(isVector(R))
        self.radius2 = np.abs(r2)
        self.R2 = normalise(isVector(R))

        self.thickness = thickness
        self.ref_index = ref_index
        self.w1 = self.radius1 - np.sqrt(self.radius1**2 - self.height**2)
        self.w2 = self.radius2 - np.sqrt(self.radius2**2 - self.height**2)
        self.w = np.mean(self.w1 + self.w2)
        self.d = self.w1+self.w2+self.thickness

        # Build the composite
        self.add(SphericalRefraction(self.centre - (self.radius1 - self.w1 - 0.5 * self.thickness) * self.R1,
                                     self.radius1 * self.R1, height=self.height, n1=air, n2=ref_index))
        if self.thickness > 0:
            self.add(CylindricalRefraction(self.centre, self.R1 * self.thickness * 0.5,
                                           self.height, self.ref_index,air, ))

        self.add(SphericalRefraction(self.centre - (self.radius2 - self.w2 - 0.5 * self.thickness) * self.R2,
                                     self.radius2 * self.R2, height=self.height, n1=air, n2=ref_index))


        # Analytic calculation of focal point
        self.f = lambda wl: lensmakers(self.radius1, -1. * self.radius2, self.ref_index(wl), self.thickness + self.w1+self.w2)
        self.bfd = lambda wl : self.f(wl) * (1 - (self.ref_index(wl)-1)*self.d/(self.ref_index(wl)*np.abs(self.radius1)) )


class BiConcave(CompositeObject):

    """BiConcave lens - Very untested"""

    def __init__(self, centre, height, R1, R2, thickness, ref_index):
        super(BiConcave, self).__init__(centre)
        self.height = height
        self.radius1 = np.linalg.norm(R1)
        self.R1 = normalise(isVector(R1))
        self.radius2 = np.linalg.norm(R2)
        self.R2 = normalise(isVector(R2))
        self.ref_index = ref_index
        self.w1 = self.radius1 - np.sqrt(self.radius1**2 - self.height**2)
        self.w2 = self.radius2 - np.sqrt(self.radius2**2 - self.height**2)
        self.thickness = thickness + self.w1 + self.w2
        self.w = np.mean(self.w1 + self.w2)

        self.add(SphericalRefraction(self.centre + (self.radius1 - 0.5 * self.w + 0.5 * self.thickness) * self.R1,
                                     -self.radius1 * self.R1, height=self.height, n1=air, n2=ref_index))

        self.add(CylindricalRefraction(self.centre + self.R1 * 0.5 * (self.w1 - self.w2), self.R1 * self.thickness * 0.5,
                                       self.height, air, self.ref_index))

        self.add(SphericalRefraction(self.centre + (self.radius2 - 0.5 * self.w + 0.5 * self.thickness) * self.R2,
                                     -self.radius2 * self.R2, height=self.height, n1=air, n2=ref_index))

    def lensmaker(self, wl=500):
        return lensmakers(-self.radius1, 1. * self.radius2, self.ref_index(wl), self.thickness + self.w1 + self.w2) + self.centre[0]
