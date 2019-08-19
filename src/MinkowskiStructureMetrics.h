#ifndef MINKOWSKI_STRUCTURE_H
#define MINKOWSKI_STRUCTURE_H

#include <voro++.hh> // Voro++ handles the Voronoi construction part
#include <complex>
#include <iostream>
#include <malloc.h>
// #include <time.h> // Just to seed RNG with local time

/******************************************************************************
  These files define an interface for calculating the Minkowski structure
	metrics as defined in Mickel et al. (2013)'s paper "Shortcomings of the
	bond orientational order parameters for the analysis of disordered
	particulate matter" (DOI: 10.1063/1.4774084). These structure metrics
	characterize the structure of a set of points in three dimensions in
	terms of a number of order parameters {q'} that are closely related to
	the Steinhardt order parameters (DOI: 10.1103/PhysRevB.28.784).
******************************************************************************/

// NOTE
// For even l spherical harmonics we don't need a bond _direction_, so can optimize those?
// If particle is anisotropic, should we rotate angles into particle frame?
//   => Probably not, since spherical harmonics are rotationally invariant.


// void msm_prepare(unsigned int n_positions, double **positions, double box[9]) // calculates neighbours
// void calc_ql(unsigned int p_idx, unsigned int l) // calculates q_l for particle p
// void calc_wl(unsigned int p_idx, unsigned int l) // calculates w_l for particle p
// void calc_all_ql(unsigned int p_idx, unsigned int l, float *q) // calculates q_l for all particles
// void calc_all_wl(unsigned int p_idx, unsigned int l, float *w) // calculates w_l for all particles




namespace MSM {


  struct Voro_Nbs {
		std::vector<std::vector<int>> indices;
		std::vector<std::vector<double>> face_areas;
		std::vector<double> cell_areas;
    std::vector<std::vector<double>> thetas;
    std::vector<std::vector<double>> phis;
  };


	// Calculates the associated Legendre polynomial P_l^m(x)
	// RvD: Checked outputs, matches GNU scientific library with rel. error ~1E-8.
	float LegendrePlm_m_gtr_0(int l, int m, double x);

	// Calculates ln(x!) so we can calculate (l-m)!/(l+m)!
	// RvD: Checked outputs, matches GNU scientific library with rel. error ~1E-16.
	double gammln(double xx);

	// Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
	double spherical_harmonic_factor(unsigned int l, int m);

	// Calculates a spherical harmonic Y_l^m(z,phi)
	std::complex<double> spherical_harmonic(int l, int m, double theta, double phi);

// input: set of coordinates in a box of which we want the Minkowski structure metrics
// output: array of requested MSM

	// C form
	void msmc(size_t n_positions, double **positions, double box[9], double (*q)[6]);

	// C++ form
	void msm(
		const std::vector<std::vector<double>>& positions,
		const std::vector<double>& box, // column-major, upper triangular (a,b,c,0,e,f,0,0,g)
		std::vector<std::vector<double>>& q,
		std::vector<std::vector<double>>& w
	);




	// Class that can calculate return specific Minkowski structure metrics (ql's).
	// Should store some values (e.g. spherical harmonic factors) for efficient reuse.
	// Interface should have various levels for certain easy specific cases (e.g.
	// spheres, which can use lab frame instead of particle frame angles).
	// class MSMCalculator {
	// private:
	// public:
	// 	MSMCalculator();
	// 	virtual ~MSMCalculator();
	// };

}





#endif
