#!/usr/bin/env python

import numpy as np


def isVector(x):
    if isinstance(x, (list, tuple, np.ndarray)) and len(x) != 3:
        raise Exception("Invalid 3d vector")
    return(np.array(x))
    
def normalise(n):
    """docstring for normalise"""
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
