A small C++ code that can read in configurations of particle coordinates and calculate the corresponding Minkowski structure metrics (MSM). 

These metrics, defined in Ref. [1][1], describe the symmetry of where a particle's neighbours are located around it. They are a variant of the Steinhardt bond order parameters [2][2] in which the neighbours are defined through a Voronoi construction. This provides a parameter-free way of defining the neighbours. Additionally, the contribution of each bond is weighted by the area of the shared Voronoi facet, which ensures that nearby bonds contribute more, while distant bonds contribute less. Furthermore, this weighting makes sure that the MSM change smoothly with bond distance, instead of discontinuously like with a fixed cutoff or a fixed number of neighbours.

Currently, the code calculates:


* q<sub>lm</sub>

![qlm](https://latex.codecogs.com/png.latex?q_{lm}(i)&space;=&space;\sum_{f\in\mathcal{F}(i)}&space;\frac{A(f)}{A}&space;Y_{lm}(\mathbf{r}_{ij}))

* averaged q<sub>lm</sub>

![avqlm](https://latex.codecogs.com/png.latex?\bar{q}_{lm}(i)&space;=&space;\frac{1}{\tilde{N}_b(i)}&space;\sum_{k=0}^{\tilde{N}_b(i)}&space;q_{lm}(k))

* q<sub>l</sub>

![ql](https://latex.codecogs.com/png.latex?q_{lm}(i)&space;=&space;\frac{1}{N_b(i)}\sum_{j=0}^{N_b(i)}&space;Y_{lm}(\mathbf{r}_{ij}))

* w<sub>l</sub>

![wl](https://latex.codecogs.com/png.latex?w_l(i)&space;=&space;\sum_{m_1&plus;m_2&plus;m_3=0}&space;\begin{pmatrix}&space;l&space;&&space;l&space;&&space;l&space;\\&space;m_1&space;&&space;m_2&space;&&space;m_3&space;\end{pmatrix}&space;q_{lm_1}(i)&space;q_{lm_2}(i)&space;q_{lm_3}(i))

* averaged q<sub>l</sub>

![avql](https://latex.codecogs.com/png.latex?\bar{q}_l(i)=\sqrt{\frac{4\pi}{2l&plus;1}\sum_{m=-l}^l&space;\left|&space;\bar{q}_{lm}(i)&space;\right|^2})

* averaged w<sub>l</sub>

![avwl](https://latex.codecogs.com/png.latex?\bar{w}_l(i)&space;=&space;\sum_{m_1&plus;m_2&plus;m_3=0}&space;\begin{pmatrix}&space;l&space;&&space;l&space;&&space;l&space;\\&space;m_1&space;&&space;m_2&space;&&space;m_3&space;\end{pmatrix}&space;\bar{q}_{lm_1}(i)&space;\bar{q}_{lm_2}(i)&space;\bar{q}_{lm_3}(i))

This code currently supports reading in both text coordinate files (with extension ".dat") and the Glotzerlab general simulation data format (".gsd"/".GSD") [3][3]. The format of the text files should be:

	N
	a1x a2x a3x a1y a2y a3y a1z a2z a3z
	p1x p1y p1z
	...
	pNx pNy pNz

where `N` is the number of particles, `a1`, `a2` and `a3` are the lattice vectors and `p1`-`pN` are particle positions. Other data may follow the particle coordinates on the same line, but it will be ignored.

GSD files should be supplied as-is and should require no tweaking.

[1]: [Mickel et al., J. Chem. Phys. 2013](scitation.aip.org/content/aip/journal/jcp/138/4/10.1063/1.4774084)
[2]: [Steinhardt et al., Phys. Rev. B 1983](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784)
[3]: [https://gsd.readthedocs.io/]
