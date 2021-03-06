# What is this?

A small C++ code that can read in configurations of 3D particle coordinates and calculate the corresponding Minkowski structure metrics (MSM).

These metrics, defined in Ref. [[1]](#Mickel2013), describe the symmetry of where a particle's neighbours are located around it. They are a variant of the Steinhardt bond order parameters [[2]](#Steinhardt1983) in which the neighbours are defined through a Voronoi construction. This provides a parameter-free way of defining the neighbours. Additionally, the contribution of each bond is weighted by the area of the shared Voronoi facet, which ensures that nearby bonds contribute more, while distant bonds contribute less. Furthermore, this weighting makes sure that the MSM change smoothly with bond distance, instead of discontinuously like with a fixed cutoff or a fixed number of neighbours. 

This particular code takes special care to make sure the method works correctly for triclinic and very small periodic boxes, something I found lacking in other implementations.

Currently, the code calculates:

* q<sub>lm</sub>
* averaged q<sub>lm</sub>
* q<sub>l</sub>
* w<sub>l</sub>
* averaged q<sub>l</sub>
* averaged w<sub>l</sub>

# Compilation

Before you can compile this code, you will need:
* [CMake](https://cmake.org/) to build the compilation instructions
* The C++ [Eigen](http://eigen.tuxfamily.org/) library
* The C++ [Voro++](http://math.lbl.gov/voro++/) library

Once you have all these, the code should be compilable.

#### Linux

Navigate to `BOOP_per_snapshot/bin` and type `cmake ..; make`. The code should build.

#### Windows / Mac

TODO. Haven't figured this out yet.

# Usage

This code currently supports reading in both text coordinate files (with extension ".dat") and the Glotzerlab general simulation data format (".gsd"/".GSD") [[3]](#GSD). The format of the text files should be:

	N
	a1x a2x a3x a1y a2y a3y a1z a2z a3z
	p1x p1y p1z
	...
	pNx pNy pNz

where `N` is the number of particles, `a1`, `a2` and `a3` are the lattice vectors and `p1`-`pN` are particle positions. Other data may follow the particle coordinates on the same line, but it will be ignored.

GSD files should be supplied as-is and should require no tweaking.

<a name="Mickel2013">1</a>:  [Mickel et al., J. Chem. Phys. (2013)](scitation.aip.org/content/aip/journal/jcp/138/4/10.1063/1.4774084)

<a name="Steinhardt1983">2</a>:  [Steinhardt et al., Phys. Rev. B (1983)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784)

<a name="GSD">3</a>:  https://gsd.readthedocs.io/
