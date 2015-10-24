#!/usr/bin/env python
from .utils import isVector, normalise,wavelengthToHex
import numpy as np
from .optical_elements import *

air = lambda l : 1.

def lensmakers(r1, r2, n, d):
    f = 1./( (n-1)*( 1./r1 - 1./r2 + (n-1)*d/(n*r1*r2)) )
    return f

class PlanoConvex(CompositeObject):
    """docstring for PlanoConvex"""
    def __init__(self, centre, height, R, thickness, ref_index):
        super(PlanoConvex, self).__init__(centre)
        self.height = height
        self.radius = np.linalg.norm(R)
        self.R = normalise(isVector(R))
        self.thickness = thickness
        self.ref_index = ref_index
        self.w = self.radius - np.sqrt(self.radius**2 - self.height**2)
        
        self.add(SphericalRefraction(self.centre - (self.radius - self.w-0.5*self.thickness)*self.R,
                                     self.radius*self.R, height=self.height, n1=air, n2=ref_index))
        
        self.add(CylindricalRefraction(self.centre,self.R*self.thickness*0.5,
                                 self.height, air, self.ref_index))
        
        a = self.height*np.array([0,0,1])
        b = self.height*normalise(np.cross(a,self.R))
        self.add(PlaneRefraction(self.centre-self.R*0.5*self.thickness,a,b, n1=air, n2=ref_index ))
        
        self.f = lambda wl : lensmakers(float('inf'), -1.*self.radius, self.ref_index(wl), self.thickness+self.w)+self.centre[0]+0.5*self.thickness+self.w


class PlanoConcave(CompositeObject):
    """docstring for PlanoConcave"""
    def __init__(self, centre, height, R, thickness, ref_index):
        super(PlanoConcave, self).__init__(centre)
        self.height = height
        self.radius = np.linalg.norm(R)
        self.R = normalise(isVector(R))
        self.ref_index = ref_index
        self.w = self.radius - np.sqrt(self.radius**2 - self.height**2)
        self.thickness = thickness + self.w
        
        self.add(SphericalRefraction(self.centre + (self.radius - self.w + 0.5*self.thickness)*self.R,
                                     -self.radius*self.R, height=self.height, n1=ref_index, n2=air))
        
        self.add(CylindricalRefraction(self.centre,self.R*self.thickness*0.5,
                                 self.height, air, self.ref_index))
        
        a = self.height*np.array([0,0,1])
        b = self.height*normalise(np.cross(a,self.R))
        self.add(PlaneRefraction(self.centre-self.R*0.5*self.thickness,a,b, n1=air, n2=ref_index ))
    
    def lensmaker(self, wl = 500):
        return lensmakers(float('inf'), 1.*self.radius, self.ref_index(wl), self.thickness+self.w)+self.centre[0]
                                   


class BiConvex(CompositeObject):
    """docstring for BiConvex"""
    def __init__(self, centre, height, R1, R2, thickness, ref_index):
        super(BiConvex, self).__init__(centre)
        self.height = height
        self.radius1 = np.linalg.norm(R1)
        self.R1 = normalise(isVector(R1))
        self.radius2 = np.linalg.norm(R2)
        self.R2 = normalise(isVector(R2))
        self.thickness = thickness
        self.ref_index = ref_index
        self.w1 = self.radius1 - np.sqrt(self.radius1**2 - self.height**2)
        self.w2 = self.radius2 - np.sqrt(self.radius2**2 - self.height**2)
        self.w = np.mean(self.w1+self.w2)
        
        self.add(SphericalRefraction(self.centre - (self.radius1 - 0.5*self.w-0.5*self.thickness)*self.R1,
                                     self.radius1*self.R1, height=self.height, n1=air, n2=ref_index))
        
        self.add(CylindricalRefraction(self.centre - self.R1*0.5*(self.w1-self.w2),self.R1*self.thickness*0.5,
                                 self.height, air, self.ref_index))
        
        self.add(SphericalRefraction(self.centre - (self.radius2 - 0.5*self.w - 0.5*self.thickness)*self.R2,
                                     self.radius2*self.R2, height=self.height, n1=air, n2=ref_index))
    
    def lensmaker(self, wl = 500):
        return lensmakers(self.radius1, -1.*self.radius2, self.ref_index(wl), self.thickness+self.w1+self.w2)+self.centre[0]




class BiConcave(CompositeObject):
    """docstring for BiConcave"""
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
        self.w = np.mean(self.w1+self.w2)
        
        self.add(SphericalRefraction(self.centre + (self.radius1 - 0.5*self.w+0.5*self.thickness)*self.R1,
                                     -self.radius1*self.R1, height=self.height, n1=air, n2=ref_index))
        
        self.add(CylindricalRefraction(self.centre + self.R1*0.5*(self.w1-self.w2),self.R1*self.thickness*0.5,
                                 self.height, air, self.ref_index))
        
        self.add(SphericalRefraction(self.centre + (self.radius2 - 0.5*self.w + 0.5*self.thickness)*self.R2,
                                     -self.radius2*self.R2, height=self.height, n1=air, n2=ref_index))
    
    def lensmaker(self, wl = 500):
        return lensmakers(-self.radius1, 1.*self.radius2, self.ref_index(wl), self.thickness+self.w1+self.w2)+self.centre[0]


