#ifndef MINKOWSKI_STRUCTURE_H
#define MINKOWSKI_STRUCTURE_H

#include <voro++.hh> // Voro++ handles the Voronoi construction part
#include <complex>
#include <iostream>
#include <malloc.h>
#include <cmath>
#include <unordered_map>
#include <eigen3/Eigen/Geometry>
#include <algorithm>
// #include <time.h> // Just to seed RNG with local time

/******************************************************************************
  These files define a class for calculating the Minkowski structure
	metrics as defined in Mickel et al. (2013)'s paper "Shortcomings of the
	bond orientational order parameters for the analysis of disordered
	particulate matter" (DOI: 10.1063/1.4774084). These structure metrics
	characterize the structure of a set of points in three dimensions in
	terms of a number of order parameters {q'} that are closely related to
	the Steinhardt order parameters (DOI: 10.1103/PhysRevB.28.784).
******************************************************************************/

namespace MSM {

struct neighbourBonds {
	std::vector<int> indices;
	std::vector<double> face_areas;
	float cell_area;
  std::vector<float> thetas;
  std::vector<float> phis;
};

class MinkowskiStructureCalculator {
private:
  std::vector<std::vector<float>> positions_;
  std::vector<float> box_;
  std::vector<neighbourBonds> nbs_;

  // Applies periodic boudary conditions
  Eigen::Vector3d nearest_image(
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Matrix3d& box
  );
  // Calculates the associated Legendre polynomial P_l^m(x)
  float LegendrePlm_m_gtr_0(int l, int m, float x);
  // Unit test for the LegendrePlm_m_gtr_0 function
  void verify_LegendreP(void);
  // Calculates ln((x-1)!) for integers x > 0 and tabulates them for reuse.
  float gammln_i(int i);
  // Calculate the Wigner 3j-symbols using the Racah formula
  float wigner3j(int l, int m1, int m2, int m3);
  // Unit test for the wigner3j function
  void verify_wigner3j(void);
  // Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
  float spherical_harmonic_factor(unsigned int l, int m);
  // Calculates a spherical harmonic Y_l^m(z,phi)
  std::complex<float> spherical_harmonic(int l, int m, float theta, float phi);
  // Unit check for the spherical_harmonic function
  void verify_spherical_harmonic(void);
  // Does a number of checks to ensure the Voro++ output is usable
  void verify_voro_results(void);
  // Copies the input configuration a few times if it is too small
  void enlarge_configuration(void);

public:
  MinkowskiStructureCalculator(void);
  virtual ~MinkowskiStructureCalculator(void){};

  // Loads in a configuration
  void msm_prepare(
    const std::vector<std::vector<float>>& positions,
    const std::vector<float>& box
  );
  // Calculates a single q_l for particle p
  float ql(unsigned int p, unsigned int l);
  // Calculates a single w_l for particle p
  float wl(unsigned int p, unsigned int l);
  // Fills the q and w vectors with their values for input positions, starting from l=0
  void compute(
      const std::vector<std::vector<float>>& positions,
      const std::vector<float>& box,
      std::vector<std::vector<float>>& q,
      std::vector<std::vector<float>>& w
  );
};


} // End namespace MSM





#endif
