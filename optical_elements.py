from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from .utils import isVector, normalise
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
from . import n_air

eps = np.finfo(np.float64).eps

# TODO Add kwargs to draw that pass to Polygon to allow for more complex rendering ie colors for materials
class OpticalElement(object):
    
    """Abstract base class defining the common interface that all optical elements 
        must implement and subclass.
        
        Also defines drawBench() which should be called first in all derived drawBench 
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
    def drawBench(self, ax=None):
        if ax == None:
            fig = plt.figure()
            ax = fig.add_axes([0.05,0.05, 0.9, 0.9], aspect='auto')
            ax.set_xlim(-2,2)
            ax.set_ylim(-2,2)
            

        return ax



class Cube(OpticalElement):
    """docstring for Cube"""
    def __init__(self, centre,  a, b, c, n):
        super(Cube, self).__init__(centre)
        self.a = isVector(a)
        self.b = isVector(b)
        self.c = isVector(c)
        
        self.n = n
        
        self.side_centres = [self.centre + self.c,
                              self.centre - self.c,
                              self.centre + self.b,
                              self.centre - self.b,
                              self.centre + self.a,
                              self.centre - self.a]
        
        self.sides = [PlaneInterface(self.side_centres[0],  self.a, self.b, self.n),
                      PlaneInterface(self.side_centres[1],  -self.a, self.b, self.n),
                      PlaneInterface(self.side_centres[2],  -self.a, self.c, self.n),
                      PlaneInterface(self.side_centres[3],  self.a, self.c, self.n),
                      PlaneInterface(self.side_centres[4],  self.b, self.c, self.n),
                      PlaneInterface(self.side_centres[5],  -self.b, self.c, self.n)]
    
    def drawBench(self, ax=None):
        """Draws a projection of the object on the xy plane"""
        ax = super(Cube, self).drawBench(ax)
        points = np.array([x[:2] for y in self.sides for x in y.points ])
        hull = ConvexHull(points)
        
        for simplex in hull.simplices:
         ax.plot(points[simplex,0], points[simplex,1], 'k-')

    def distance(self, ray):
        d = np.array([s.distance(ray) for s in self.sides])
        
        # Filter out negative distances but leave 0 distance
        d[np.logical_or(np.isnan(d), d < 0)] = float('inf')
        
        self.entrance_idx = np.argmin(d) 
        return d
                
    def propagate_ray(self, ray):  
        # Propagates the ray into the block
        ray = self.sides[self.entrance_idx].propagate_ray(ray)
        
        ## Decides which side it will exit through (can't be entrance side)
        ## TODO this is where TIR could be incorporated
        d = np.array([s.distance(ray) for s in self.sides])
        d[np.logical_or(np.isnan(d), d <= 0+eps)] = float('inf')
        self.exit_idx = np.argmin(d)
        ray.append(ray.p + np.min(d)*ray.k, ray.k)
        
        return ray
        
    def exitToAir(self,ray):
        """docstring for exitToAir"""
        ray = self.sides[self.exit_idx].exitToAir(ray)
        return ray
class Spherical(OpticalElement):
    """docstring for Spherical"""
    
    def __init__(self, centre, R, theta):
        super(Spherical, self).__init__(centre)
        
        # scalar radius
        self.r = np.linalg.norm(R)
        # Unit vector defining radial direction
        self.R = normalise(R)    
        
        self.theta = theta
        self.t_0 = np.arctan2(self.R[1], self.R[0])

        assert np.logical_and( self.theta > 0, theta <= np.pi ), "0 < theta <= pi"
        
        if self.theta < 0.5*np.pi:
            self.height = self.r * np.sin(self.theta)
        else:
            self.height = self.r


    def drawBench(self, ax=None):
        """docstring for drawBench"""
        ax = super(Spherical, self).drawBench(ax)
        
        
        if self.R[2] == 0:
            # If the spherical isn't tilted in the z direction we can 
            # use a simplified plotting algorithm since we know the boundary
                    
            T = np.linspace(self.t_0-self.theta, self.t_0+self.theta, 100)
            points = [[self.r*np.cos(t)+self.centre[0], self.r*np.sin(t)+self.centre[1]] for t in T]
            poly = Polygon(points, False, facecolor = 'none')
            ax.add_patch(poly)
        
        # TODO Make it so that the spherical can have a z comp
        
        else:
            raise NotImplementedError()
        #    # Otherwise use a more complex algorithm to find the convex hull
        #    # of the 2d projection of the shape
        #    T = np.linspace(self.t_0-self.theta, self.t_0+self.theta, 100)
        #    P = np.linspace(0, 2*np.pi, 100)
        #    
        #    points = np.array([[self.r*np.sin(t1)*np.cos(t2)+self.centre[0], 
        #                        self.r*np.sin(t1)*np.sin(t2)+self.centre[1]] for t in T for p in P])
        #    hull = ConvexHull(points)
        #    
        #    for simplex in hull.simplices:
        #        ax.plot(points[simplex,0], points[simplex,1], 'k-')
        #    #ax.scatter(points[:,0], points[:,1],points[:,2])
            
    def plot3D(self):
        """docstring for plot3D"""
        raise NotImplementedError()
        
        #fig = plt.figure()
        #ax = fig.add_axes([0.05,0.05, 0.9, 0.9], projection='3d')
        #ax.set_xlim(-1.5,1.5)
        #ax.set_ylim(-1.5,1.5)
        #ax.set_zlim(-1.5,1.5)
        #
        #T = np.linspace(self.t_0-self.theta, self.t_0+self.theta, 100)
        #P = np.linspace(self.p_0-self.phi, self.p_0+self.phi, 100)
        #
        #points = np.array([[self.r*np.sin(t)*np.cos(p)+self.centre[0], 
        #                    self.r*np.sin(t)*np.sin(p)+self.centre[1],
        #                    self.r*np.cos(t)+self.centre[2]] for t in T for p in P])
        #
        #ax.scatter(points[:,0], points[:,1],points[:,2])

    def distance(self, ray):
        r = ray.p - self.centre
        # TODO Fix this to mach the new list return type
        # If the sqrt is imaginary it never intersects
        sq = np.dot(r, ray.k)**2 - (np.linalg.norm(r)**2 - self.r**2)
        if sq >= 0 :
            l1 = -np.dot(r, ray.k) + np.sqrt(sq)
            l2 = -np.dot(r, ray.k) - np.sqrt(sq)
        else:
            return [float('inf')]
        
        # Try the closest intersection first
        d = sorted([l1,l2]) 
        
        # Note using machine eps to to sqrts 
        # this was a bug
        if d[0] > eps and self.isIntersect(ray.p + d[0]*ray.k - self.centre):
            return [d[0]] 
        # Now try the further away one                
        elif d[1] > eps and self.isIntersect(ray.p + d[1]*ray.k - self.centre):
            return [d[1]]        
        else: 
            return [float('inf')]
    
    def isIntersect(self,p):
        """docstring for isIntersect"""
        
        ang = np.arccos(np.dot(p,self.R)/np.linalg.norm(p))
        z_cond = np.abs(p[2]) <= self.height
        return np.logical_and(z_cond, np.abs(ang) <= self.theta)



class SphericalWall(Spherical):
    """A spherical with total absorption"""
    def __init__(self, centre,  R, theta):
        super(SphericalWall, self).__init__( centre, R, theta)
    
    def propagate_ray(self, ray):
        """Implements propagate_ray for Wall by terminating ray"""
        d = self.distance(ray)        
        ray.append(ray.p + d*ray.k, ray.k)
        ray.isTerminated = True
        return ray       

class SphericalMirror(Spherical):
    """A totally reflective spherical (both sides)"""
    def __init__(self, centre,  R, theta):
        super(SphericalMirror, self).__init__( centre, R, theta)
    
    def propagate_ray(self, ray):
        """Implements propagate_ray for Wall by terminating ray"""

        d = self.distance(ray) 
        p = ray.p + d*ray.k  
        
        # Normal at the point of intersection
        n = normalise(p - self.centre)  
        
        # Orthogonal to both k and n
        m = normalise(np.cross(n, np.cross(ray.k, n)))
        
        
        # New wavevector
        #k_prime = np.dot(ray.k, m)*m - np.dot(ray.k, n)*n
         
        k_prime = ray.k - 2*np.dot(ray.k, n)*n        
           
        ray.append(p, k_prime)
        ray.isTerminated = False
        return ray    




class Plane(OpticalElement):
    """Plane optical element abstract class, used for OpticalBench 
        boundaries, mirrors etc
       
        pos - Vector to the center of the plane
        a,b - Vectors describing the distance from the center to each edge,
              *must be orthogonal*. Side lengths are then 2a, 2b"""
    
    def __init__(self, centre,  a, b):
        super(Plane, self).__init__(centre)
        
        assert np.dot(a, b) == 0, "Vectors must be orthogonal"
        
        self.a = isVector(a)
        self.b  = isVector(b)
        
        n = np.cross(a,b)
        self.normal = normalise(np.cross(a,b))
        
        # The four corners of the plane
        self.points = [self.centre + a + b,
                       self.centre - a + b,
                       self.centre - a - b,
                       self.centre + a - b ]
        
        
    def drawBench(self, ax=None):
        """Draws a projection of the object on the xy plane"""
        ax = super(Plane, self).drawBench(ax)
        poly = Polygon([x[:2] for x in self.points], True, facecolor = 'none')
        ax.add_patch(poly)
        
    
    def distance(self, ray):
        """Distance from ray to self along the ray's path"""    
        
        r = self.centre - ray.p
        d = np.dot(self.normal, r)/np.dot(ray.k, self.normal) #if np.dot(self.normal, r) != 0 else float('inf')
        
        p = ray.p + d*ray.k
        d = d if self.isIntersect(p) else float('inf')
        d = float('inf') if np.isnan(d) or d < 0 else d
        return [d]
        
    def isIntersect(self, p):
        """Decides if p is in the plane or not by projecting the intersection point
            on the infinite plane onto the a,b axis than span the finite plane
            and seeing if either component is >1"""
        alpha = np.dot(self.a, (p-self.centre))/np.linalg.norm(self.a)**2
        beta = np.dot(self.b, (p-self.centre))/np.linalg.norm(self.b)**2
        
        return np.logical_and(np.abs(alpha) <= 1, np.abs(beta) <= 1)
        
class Mirror(Plane):
    """A totally reflective plane (both sides)"""
    def __init__(self, centre,  a, b):
        super(Mirror, self).__init__(centre,  a, b)
    
    def propagate_ray(self, ray):
        """Implements propagate_ray for Mirror by reflecting the wavevector of the ray
            about the normal of the plane.
            """
        k_prime = ray.k - 2*np.dot(ray.k, self.normal)*self.normal 
        d = self.distance(ray)
        ray.append(ray.p + d*ray.k, k_prime)        
        return ray
        
class Wall(Plane):
    """A plane with total absorption"""
    def __init__(self, centre,  a, b):
        super(Wall, self).__init__(centre,  a, b)
    
    def propagate_ray(self, ray):
        """Implements propagate_ray for Wall by terminating ray"""
        d = self.distance(ray)        
        ray.append(ray.p + d*ray.k, ray.k)
        ray.isTerminated = True
        return ray
        

class Screen(Plane):
    """An output screen, records all incident beams as pixels and renders them"""
    def __init__(self, centre,  a, b):
        super(Screen, self).__init__(centre,  a, b)
        
        self.fig = plt.figure()
        self.ax = self.fig.add_axes([0.05,0.05, 0.9, 0.9], aspect='auto')
        A = np.linalg.norm(self.a)
        B = np.linalg.norm(self.b)
        
        self.ax.set_xlim(-1*A, A)
        self.ax.set_ylim(-1*B, B)
        self.ax.set_aspect('equal')
        
    
    def propagate_ray(self, ray):
        """Implements propagate_ray for Screen, projects the ray's intersection with
            the screen onto plane orthogonal coordinates and draws it as a point on the
            associated mpl figure"""    
        d = self.distance(ray)        
        ray.append(ray.p + d*ray.k, ray.k)
        ray.isTerminated = True
        
        alpha = np.dot(self.a, (ray.p-self.centre))/np.linalg.norm(self.a)
        beta = np.dot(self.b, (ray.p-self.centre))/np.linalg.norm(self.b)
        
        self.ax.scatter(alpha, beta)
        return ray
        



class PlaneInterface(Plane):
    """docstring for PlaneInterface"""
    def __init__(self, centre,  a, b,n2):
        super(PlaneInterface, self).__init__(centre,  a, b)
        self.n2 = n2
     
    def propagate_ray(self, ray):     
        """Implements propagate_ray for Wall by terminating ray"""     
        d = self.distance(ray)     
             
        r = ray.n/self.n2    
        c = -np.dot(self.normal, ray.k) 
        if 1-(r**2)*(1-c**2) < 0:
            print("TOTAL INTERNAL REFLECTION")
        k_prime = r*ray.k + (r*c - np.sqrt(1-(r**2)*(1-c**2)))*self.normal     
        ray.append(ray.p + d*ray.k, k_prime)     
        ray.isTerminated = False  
        ray.n = self.n2   
        return ray 
    
    def exitToAir(self, ray):     

        d = self.distance(ray)     
             
        r = ray.n/n_air   
        c = -np.dot(-1*self.normal, ray.k) 
        if 1-(r**2)*(1-c**2) < 0:
            print("TOTAL INTERNAL REFLECTION")
        k_prime = r*ray.k + (r*c - np.sqrt(1-(r**2)*(1-c**2)))*-1*self.normal     
        ray.append(ray.p + d*ray.k, k_prime)     
        ray.isTerminated = False  
        ray.n = n_air  
        return ray        
        

