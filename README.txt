
Raytracer-0.0.1

Getting Started:

All simulations essentially are constructed by creating an OpticalBench object,
 adding Sources, OpticalElements and Screens to it, then calling Render(), then
 draw()

The file examples.py contains some simple examples, this should simply run as

	python examples.py
	
Example:

I recommend you import the module using

>>>import raytracer
>>>import raytracer.optical_elements as oe
>>>import raytracer.sources as src
>>>import raytracer.lenses as lens
>>>import raytracer.utils as utils

And also define the standard refractive indices using Cauchy's equation

>>>glass = lambda l : 1.5046 + 4200./l**2
>>>air = lambda l : 1.                                
>>>water =  lambda l : 1.319 + 6878./l**2 
  
Create an Optical bench with length, width and height all 2*50 = 100mm
>>>ob = raytracer.OpticalBench(50, 50, 50, verbose=True)


Add a ray, starting a coordinates (20,50,50), travelling in direction [1,0,0]
 with wavelength 450nm

>>>ob.add(src.SingleRay([20,50, 50], [1,0,0], 450))

Add a lens, centred at (50,50,50), the centre of the OpticalBench.
>>>ob.add(lens.BiConvex([50,50,50], height=22, r1=50, r2=50, R=[1,0,0],
thickness=1., ref_index=glass))


Now call render, this populates the vertices attributes of all of the rays in
the OpticalBench.

>>>ob.Render()

Finally call draw.

>>>ob.draw()
>>>plt.show()




Potential GOTCHAs:

- Almost all lengths a specified as vectors spanning half the distance.
	ie a,b for a plane point from the centre to each edge, so the total plane is 2a*2b.
	
- Planoconvex/concave lenses consist of a spherical and a plane. The plane
 sticks outside of the curved surface at the corners (it's a square enclosing a
  circle), this means a ray that catches this region will trigger an exception
   because it's trying to refract from glass to air, having not first refracted
   from air to glass


- Concave lenses are not well tested.

- Getting the focal point correct can be tricky



Unit tests:

test.py contains some unit tests based on the unittest framework.
These need to be run by executing 

python -m unittest discover

from inside the main raytracer root directory.

or by 

nosetests raytracer/test.py

If you don't have unittest/nose installed the tests won't run. But they all
passed before I made my submission to Blackboard.

The test coverage is 78%, it's mostly exception paths that aren't fully tested.