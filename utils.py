#!/usr/bin/env python

import numpy as np
                         


def isVector(x):
    if isinstance(x, (list, tuple, np.ndarray)) and len(x) != 3:
        raise Exception("Invalid 3d vector")
    return(np.array(x))
    
def normalise(n):
    """docstring for normalise"""
    if np.linalg.norm(n) == 0 :
        raise Exception("Attempt to normalise 0 vector")
    return np.array(n) / np.linalg.norm(n)
    
def rtpairs(R,N):
    '''Generates n_i points uniformly distributed about a circle of 
    radius r_i for r_i in R and n_i in N'''
    
    assert len(R) == len(N), "Must supply same number of points and radii!"
    for r,N in zip(R,N):
        theta = 0
        for i in range(N):
            theta += 2*np.pi/N
            yield (r, theta)
            
def unifDisk(n, rmax, m):
    ''' Generates points uniformly distributed over a disk,
        n - number of rings
        rmax - maximum radius 
        m - point desity
        '''
    R = np.linspace(0, rmax, n)
    T = [i*m for i in range(n)]
    return rtpairs(R, T)

def rotPlane(a,b,t):
    return ([a*np.cos(t), a*np.sin(t),0], [0,0,b])

def rotCube(a,b,c,t):
    return ([a*np.cos(t), a*np.sin(t),0], [b*np.cos(t), -b*np.sin(t),0],  [0,0,c])

def wavelengthToHex(wl):
    if (wl >= 380 and wl < 440):
        R = -1 * (wl - 440.) / (440. - 380.)
        G = 0
        B = 1
    elif (wl >= 440. and wl < 490.):
        R = 0
        G = (wl - 440.) / (490 - 440.)
        B = 1  
    elif (wl >= 490 and wl < 510):
        R = 0
        G = 1
        B = -1 * (wl - 510.) / (510. - 490.)
    elif (wl >= 510 and wl < 580):
        R = (wl - 510.) / (580. - 510.)
        G = 1
        B = 0
    elif (wl >= 580. and wl < 645.):
        R = 1
        G = -1 * (wl - 645.) / (645. - 580.)
        B = 0.0
    elif (wl >= 645 and wl <= 780):
        R = 1
        G = 0
        B = 0
    else:
        R = 0
        G = 0
        B = 0
       
    if (wl > 780 or wl < 380):
        alpha = 0
    elif (wl > 700):
        alpha = (780. - wl) / (780. - 700.);
    elif (wl < 420):
        alpha = (wl - 380.) / (420. - 380.);
    else:
        alpha = 1
    return (R,G,B)
        
def validDistance(d):
    """Masks potential invalid distances with +inf"""
    if isinstance(d, (list, np.ndarray)):
        d = np.array(d)
        mask = np.logical_or(np.logical_or(np.isnan(d), np.isclose(d, 0.)), d < 0)
        d[mask] = float('inf')
        return d
    else:
        if np.isnan(d) or np.isclose(d, 0.) or d < 0: 
            return float('inf')
        else:
            return d