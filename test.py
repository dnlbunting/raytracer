#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose
import unittest
from scipy.optimize import minimize
from functools import partial
from .optical_elements import *
from .sources import Source, SingleRay, CollimatedBeam, Ray
from . import OpticalBench
from .lenses import *
from .utils import normalise
from .interactions import ReflectionMixin, RefractionMixin

eps = np.finfo(np.float64).eps
water = lambda l: 1.319 + 6878. / l**2
air = lambda l: 1.
glass = lambda l: 1.5046 + 4200. / l**2


class ReflectionTest(unittest.TestCase):

    def setUp(self):
        self.reflector = ReflectionMixin()
        self.reflector.distance = lambda x: 0

    def test_random(self):
        for i in range(10):
            k = normalise(np.random.uniform(0, 1, 3))
            p = np.random.uniform(0, 1, 3)
            ray = Ray(p, k, 500)
            n = normalise(np.random.uniform(0, 1, 3))

            ray = self.reflector._reflect(ray, n)
            self.assertAlmostEqual(np.dot(ray.k, n), -np.dot(k, n))


class RefractionTest(unittest.TestCase):

    def setUp(self):
        self.refractor = RefractionMixin()
        self.refractor.distance = lambda x: 0

    def test_random(self):
        for i in range(10):
            k = normalise(np.random.uniform(-1, 1, 3))
            p = np.random.uniform(0, 1, 3)

            normal = normalise(np.random.uniform(-1, 1, 3))

            self.refractor.n1 = lambda l: 1.2
            self.refractor.n2 = lambda l: 1.5

            curr_n = self.refractor.n1(500) if np.dot(k, normal) < 0 else self.refractor.n2(500)
            ray = Ray(p, k, 500)
            ray.n = curr_n

            incident = np.arccos(-np.dot(normal, k))
            ray = self.refractor._refract(ray, normal)
            refracted = np.arccos(-np.dot(normal, ray.k))

            # Snells law as ratio of sines
            self.assertAlmostEqual(np.sin(incident) / np.sin(refracted), ray.n / curr_n)


def spotRadius(f, init):
    ob = init(f)
    screen = ob.screen_list[0]
    r = np.sqrt(np.mean([x[0]**2+x[1]**2 for x in screen.pixels]))
    return(r)    
class LensTest(unittest.TestCase):

        
    def test_planoconvex(self):    
        def planoconvex(f, r=50., d=1., t=5.):
            ob = OpticalBench(150, 50, 50, verbose=False)
            ob.add(CollimatedBeam([1,50,50], [1,0,0], d, 5, 588))   
            ob.add(PlanoConvex([50, 50, 50], 40, [r,0,0], t, glass))
            ob.add(Screen([f,50, 50], [0,0,100], [0,100,0]))
            ob.Render()
            return(ob)
        
        init = partial(planoconvex)        
        opt = minimize(spotRadius, [200], bounds = [(55,399)], args=(init))
        f_optimal = opt.x
            
        ob = init(f_optimal)
        l = ob.element_list[0]
        
        self.assertAlmostEqual(f_optimal, l.bfd2(588)+l.centre[0]+0.5*l.thickness+l.w, 1)  
        
    def test_biconvex(self):    
        
        
        def biconvex(f, r1=50., r2=50.,t=5., d=1. ):
            ob = OpticalBench(150, 50, 50, verbose=False)
            ob.add(CollimatedBeam([1,50,50], [1,0,0], d, 5, 588))        
            ob.add(BiConvex([50, 50, 50], 40, r1,r2, [1,0,0], t, glass))
            ob.add(Screen([f,50, 50], [0,0,100], [0,100,0]))
            ob.Render()
            return(ob)
        init = partial(biconvex, r1=50., r2=75.)
        
        opt = minimize(spotRadius, [200], bounds = [(55,399)], args=(init))
        f_optimal = opt.x
            
        ob = init(f_optimal)
        l = ob.element_list[0]
        
        self.assertAlmostEqual(f_optimal, l.bfd(588)+l.centre[0]+0.5*l.thickness+l.w2, 1)
        

# Plane Tests


class PlaneRefractionTest(unittest.TestCase):

    def setUp(self):
        self.ob = OpticalBench(0.5, 0.5, .5)
        self.ob. add(PlaneRefraction([0.5, 0.5, .5], [0, 0, 0.1], [0, 0.1, 0],  air, water), )

    def test_perpendicular(self):
        Y = np.linspace(0, 1, 20)
        for y in Y:
            self.ob. add(SingleRay([0.1, y, 0.5], [1, 0, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, y in zip(rays, Y):
            if np.abs(y - 0.5) > 0.1:
                # Above/below plane, should hit far side
                assert_allclose(r.k, [1, 0, 0])
                assert_allclose(r.p, [1, y, 0.5])
            else:
                # Hit plane, should go straight through
                assert_allclose(r.k, [1, 0, 0])
                assert_allclose(r.p, [1, y, 0.5])

    def test_angles(self):

        T = np.linspace(-1, 1, 20)
        for t in T:
            self.ob. add(SingleRay([0.4, 0.5, 0.5], [1, t, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, t in zip(rays, T):
            incident = np.arccos(np.dot(self.ob.element_list[0].normal, normalise([1, t, 0])))
            refracted = np.arccos(np.dot(self.ob.element_list[0].normal, r.k))
            # Snells law as ratio of sines
            self.assertAlmostEqual(np.sin(incident) / np.sin(refracted), water(r.wavelength))


class SphericalNormalTest(unittest.TestCase):

    def test_spherical_norma(self):
        pass
        #ob = raytracer.OpticalBench(0.5,0.5,.5)
        #s = oe.SphericalMirror([0.5,0.5,.5], [-0.1,0,0], np.pi)
        # ob.add(s)
        #
        #Y = np.linspace(0,1, 20)
        #ax =ob.draw()
        #
        # for y in Y:
        #    a = [0.1,y,0.5]
        #    b = [y,0.1,0.5]
        #    c = [0.9,y,0.5]
        #    d = [y,0.9,0.5]
        #
        #    n = normalise(a - s.centre)
        #    ax.plot([n[0]+s.centre[0], s.centre[0]], [n[1]+s.centre[1], s.centre[1]])
        #
        #    n = normalise(b - s.centre)
        #    ax.plot([n[0]+s.centre[0], s.centre[0]], [n[1]+s.centre[1], s.centre[1]])
        #
        #    n = normalise(c - s.centre)
        #    ax.plot([n[0]+s.centre[0], s.centre[0]], [n[1]+s.centre[1], s.centre[1]])
        #
        #    n = normalise(d - s.centre)
        #    ax.plot([n[0]+s.centre[0], s.centre[0]], [n[1]+s.centre[1], s.centre[1]])


class SphericalReflectionTest(unittest.TestCase):

    def setUp(self):
        self.ob = OpticalBench(0.5, 0.5, .5)
        self.ob. add(SphericalMirror([1, 0.5, .5], [-1, 0, 0], 1))

    def test_perpendicular(self):
        Y = np.linspace(0, 1, 20)
        for y in Y:
            self.ob. add(SingleRay([0.1, y, 0.5], [1, 0, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, y in zip(rays, Y):
            if np.abs(y - 0.5) > self.ob.element_list[0].height:
                # Above/below mirror, should hit far side
                assert_allclose(r.k, [1, 0, 0])
                assert_allclose(r.p, [1, y, 0.5])
            # else:
            #    # Hit mirror, should bounce back
            #    assert_allclose(r.k, [-1,0,0])
            #    assert_allclose(r.p, [0,y,0.5])


class PlaneReflectionTest(unittest.TestCase):

    def setUp(self):
        self.ob = OpticalBench(0.5, 0.5, .5)
        self.ob. add(Mirror([0.5, 0.5, .5], [0, 0.1, 0], [0, 0, 0.1]))

    def test_perpendicular(self):
        Y = np.linspace(0, 1, 20)
        for y in Y:
            self.ob. add(SingleRay([0.1, y, 0.5], [1, 0, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, y in zip(rays, Y):
            if np.abs(y - 0.5) > 0.1:
                # Above/below mirror, should hit far side
                assert_allclose(r.k, [1, 0, 0])
                assert_allclose(r.p, [1, y, 0.5])
            else:
                # Hit mirror, should bounce back
                assert_allclose(r.k, [-1, 0, 0])
                assert_allclose(r.p, [0, y, 0.5])

    def test_angles(self):
        T = np.linspace(-1, 1, 20)
        for t in T:
            self.ob. add(SingleRay([0.4, 0.5, 0.5], [1, t, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, t in zip(rays, T):
            assert_allclose(np.arccos(np.dot(self.ob.element_list[0].normal, r.k)), (np.arccos(np.dot(self.ob.element_list[0].normal, r.k))), eps)
            assert_allclose(r.k, normalise([-1, t, 0]))

    def test_perpendicularReverse(self):
        Y = np.linspace(1, 1, 20)
        for y in Y:
            self.ob. add(SingleRay([0.9, y, 0.5], [-1, 0, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, y in zip(rays, Y):
            if np.abs(y - 0.5) > 0.1:
                # Above/below mirror, should hit far side
                assert_allclose(r.k, [-1, 0, 0])
                assert_allclose(r.p, [0, y, 0.5])
            else:
                # Hit mirror, should bounce back
                assert_allclose(r.k, [1, 0, 0])
                assert_allclose(r.p, [1, y, 0.5])

    def test_anglesReverse(self):
        T = np.linspace(-1, 1, 20)
        for t in T:
            self.ob. add(SingleRay([0.6, 0.5, 0.5], [-1, t, 0], 500))

        self.ob.Render()
        rays = [x for s in self.ob.source_list for x in s]
        for r, t in zip(rays, T):
            assert_allclose(np.arccos(np.dot(self.ob.element_list[0].normal, r.k)), (np.arccos(np.dot(self.ob.element_list[0].normal, r.k))), eps)
            assert_allclose(r.k, normalise([1, t, 0]))
