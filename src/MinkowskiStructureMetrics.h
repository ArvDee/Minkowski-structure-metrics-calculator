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
#include <memory>
#include <map>

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

class MinkowskiStructureCalculator {
private:
  // Struct for data to store before calculation of metrics
  struct particleData {
    std::vector<int> nb_indices;
    std::vector<double> thetas; // Bond angles
    std::vector<double> phis;   // Bond angles
    std::vector<double> nb_face_areas;
    double total_face_area;
  };

  std::vector<Eigen::Vector3d> positions_;
  Eigen::Matrix3d box_;
  std::vector<particleData> pData_;
  voro::container_periodic* con_;

  // These store calculated qlms as an associative map, with the pair (l,m) as
  // a key. This allows us to calculate every qlm only once per input config.
  std::map< std::pair<unsigned int, int>, std::vector<std::complex<double>> > all_qlms_;
  std::map< std::pair<unsigned int, int>, std::vector<std::complex<double>> > all_qlm_avs_;


  // Applies periodic boudary conditions
  Eigen::Vector3d nearest_image(
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Matrix3d& box
  )const;
  // Rotates the configuration so that its box matrix is upper triangular
  void rotate_box_to_uppertriangular(void);
  // Calculates the associated Legendre polynomial P_l^m(x)
  double LegendrePlm_m_gtr_0(int l, int m, double x)const;
  // Unit test for the LegendrePlm_m_gtr_0 function
  void verify_LegendreP(void)const;
  // Calculates ln((x-1)!) for integers x > 0 and tabulates them for reuse.
  double gammln_i(int i)const;
  // Calculate the Wigner 3j-symbols using the Racah formula
  double specialized_wigner3j(int l, int m1, int m2)const;
  // Unit test for the wigner3j function
  void verify_wigner3j(void)const;
  // Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
  double spherical_harmonic_factor(unsigned int l, int m)const;
  // Calculates a spherical harmonic Y_l^m(z,phi)
  std::complex<double> spherical_harmonic(int l, int m, double theta, double phi)const;
  // Unit check for the spherical_harmonic function
  void verify_spherical_harmonic(void)const;
  // Does a number of checks to ensure the Voro++ output is usable
  void verify_voro_results(void)const;
  // Checks whether the current system is too small to yield valid Voro++ output
  bool is_system_too_small(void)const;
  // Copies the input configuration a few times if it is too small
  void enlarge_system(void);

public:
  MinkowskiStructureCalculator(void);
  virtual ~MinkowskiStructureCalculator(void);

  // Loads in a configuration and prepares everything for ql / wl calculation
  void load_configuration(
    const std::vector<std::vector<double>>& positions,
    const std::vector<double>& a1,
    const std::vector<double>& a2,
    const std::vector<double>& a3
  );

  // Calculates a single weighted q_lm for one particle p
  std::complex<double> qlm(size_t p_idx, unsigned int l, int m)const;
  // Calculates a single weighted averaged q_lm_av for one particle p
  std::complex<double> qlm_av(size_t p_idx, unsigned int l, int m)const;
  // Calculates a single q_l for one particle p
  double ql(size_t p_idx, unsigned int l)const;
  // Calculates a single w_l for one particle p
  double wl(size_t p_idx, unsigned int l)const;
  // Calculates a single averaged q_l_av for one particle p
  double ql_av(size_t p_idx, unsigned int l)const;
  // Calculates a single averaged w_l_av for one particle p or
  double wl_av(size_t p_idx, unsigned int l)const;

  // Compute the qlm/ql/wl and their averages for all particles, storing
  // calculated values to avoid unnecessary re-calculations.
  const std::vector<std::complex<double>>& qlm_all(unsigned int l, int m);
  std::vector<double> ql_all(unsigned int l);
  std::vector<double> wl_all(unsigned int l);
  const std::vector<std::complex<double>>& qlm_av_all(unsigned int l, int m);
  std::vector<double> ql_av_all(unsigned int l);
  std::vector<double> wl_av_all(unsigned int l);

  // Another bond order parameter variant from DOI: 10.1039/C8CP05248D
  double ql_dot(size_t p_idx, unsigned int l);
  std::vector<double> ql_dot_all(unsigned int l);
  std::vector<double> ql_dot_av_all(unsigned int l);

  double bond_crystallinity(size_t i, size_t j, unsigned int l);
  double bond_crystallinity_lav(size_t i, size_t j, std::vector<unsigned int> all_l);

  // Number of Voronoi neighbours
  int n_voro_neighbours(size_t i);
  std::vector<int> n_voro_neighbours_all(void);
};

} // End namespace MSM

#endif
