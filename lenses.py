#!/usr/bin/env python
from .utils import isVector, normalise,wavelengthToHex
import numpy as np
from .optical_elements import *

air = lambda l : 1.                                
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
        
        self.add(SphericalRefraction(self.centre- (self.radius - self.w)*self.R,
                                     self.radius*self.R, height=self.height, n1=air, n2=ref_index))
        
        self.add(CylindricalRefraction(self.centre-self.R*self.thickness*0.5,self.R*self.thickness*0.5,
                                 self.height, air, self.ref_index)) 
        
        a = self.height*np.array([0,0,1])
        b = self.height*normalise(np.cross(a,self.R))
        self.add(PlaneRefraction(self.centre-self.R*self.thickness,a,b, n1=air, n2=ref_index ))                             

class BiConvex(CompositeObject):
    """docstring for BiConvex"""
    def __init__(self, centre, height, R, thickness, ref_index):
        super(BiConvex, self).__init__(centre)
        self.height = height
        self.radius = np.linalg.norm(R)
        self.R = normalise(isVector(R))
        self.thickness = thickness
        self.ref_index = ref_index
        self.w = self.radius - np.sqrt(self.radius**2 - self.height**2)
        
        self.add(SphericalRefraction(self.centre- (self.radius - self.w)*self.R,
                                     self.radius*self.R, height=self.height, n1=air, n2=ref_index))
    
        self.add(CylindricalRefraction(self.centre-self.R*self.thickness*0.5,self.R*self.thickness*0.5,
                                 self.height, air, self.ref_index)) 
        
        self.add(SphericalRefraction(self.centre+ (self.radius - self.w - self.thickness)*self.R,
                                     -self.radius*self.R, height=self.height, n1=air, n2=ref_index))
        




