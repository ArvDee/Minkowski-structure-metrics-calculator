#include "MinkowskiStructureMetrics.h"

/******************************************************************************
  These files define a class for calculating the Minkowski structure
  metrics as defined in Mickel et al. (2013)'s paper "Shortcomings of the
  bond orientational order parameters for the analysis of disordered
  particulate matter" (DOI: 10.1063/1.4774084). These structure metrics
  characterize the structure of a set of points in three dimensions in
  terms of a number of order parameters {q'} that are closely related to
  the Steinhardt order parameters (DOI: 10.1103/PhysRevB.28.784).
******************************************************************************/

// Empty namespace to locally define some unit test functions
namespace {
  // Verify that value is identical to ref_value within tolerance. If not, raise an error.
  void verifySimilarity(
    const double& value,
    const double& ref_value,
    const double tolerance,
    const std::string location
  ){
    double diff = value - ref_value;
    double relative_diff = diff / ref_value;
    if(fabs(relative_diff) < tolerance){
      return;
    } else { // Doing it this way also catches diff=NaN
      std::cerr << "Fatal error in Minkowski metric, location '"<<location<<"'\n";
      std::cerr << "Value, reference value, relative difference and tolerance are:\n";
      std::cerr << value <<' '<< ref_value <<' '<< relative_diff <<' '<< tolerance <<'\n';
      exit(42);
    }
  }
  // Overload function to also handle the complex values from the spherical harmonics
  void verifySimilarity(
    const std::complex<double>& value,
    const std::complex<double>& ref_value,
    const double tolerance,
    const std::string location
  ){
    std::complex<double> diff = value - ref_value;
    double relative_diff = abs(diff) / abs(ref_value);
    if(fabs(relative_diff) < tolerance){
      return;
    } else { // Doing it this way also catches diff=NaN
      std::cerr << "Fatal error in Minkowski metric, location '"<<location<<"'\n";
      std::cerr << "Value, reference value, relative difference and tolerance are:\n";
      std::cerr << value.real() <<"+"<< value.imag() <<"i ";
      std::cerr << ref_value.real() <<"+"<< ref_value.imag() <<"i ";
      std::cerr << relative_diff <<' '<< tolerance <<'\n';
      exit(42);
    }
  }

}

namespace MSM {

  // Constructor that performs some cheap unit test upon initialization
  MinkowskiStructureCalculator::MinkowskiStructureCalculator(void){
    verify_spherical_harmonic();
    verify_LegendreP();
    verify_wigner3j();
  }

  // Destructor.
  MinkowskiStructureCalculator::~MinkowskiStructureCalculator(void){
    if(con_) delete con_;
  }

  // Applies periodic boudary conditions
  Eigen::Vector3d MinkowskiStructureCalculator::nearest_image(
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Matrix3d& box
  )const{
    // Positions in relative coordinates
    Eigen::Vector3d rpos1 = box.inverse() * pos1;
    Eigen::Vector3d rpos2 = box.inverse() * pos2;
    // printf("%f %f %f\n",rpos1[0],rpos1[1],rpos1[2]);
    // printf("%f %f %f\n",rpos2[0],rpos2[1],rpos2[2]);
    // Relative non-periodic distance vector
    Eigen::Vector3d r12 = rpos2 - rpos1;
    // Absolute distance vector, to be updated in the next few lines
    Eigen::Vector3d nearvec = pos2 - pos1;
    // Check each cardinal direction for distance > box length, then correct nearvec accordingly
    for(size_t d = 0; d < 3; d++){
      if(     r12[d] > +0.5){r12[d] -= 1.0; nearvec -= box.col(d);}
      else if(r12[d] < -0.5){r12[d] += 1.0; nearvec += box.col(d);}
    }
    // printf("%f %f %f\n",r12[0],r12[1],r12[2]);
    // Return what is now the nearest image vector
    return nearvec;
  }

  // Rotates the configuration so that its box matrix is upper triangular
  void MinkowskiStructureCalculator::rotate_box_to_uppertriangular(void){
    Eigen::Vector3d a1 = box_.col(0), a2 = box_.col(1), a3 = box_.col(2);
    Eigen::Vector3d xhat, yhat, zhat;
    xhat << 1.0, 0.0, 0.0;
    yhat << 0.0, 1.0, 0.0;
    zhat << 0.0, 0.0, 1.0;

    // Find quaternion to rotate first lattice vector to x axis
    double a1mag = a1.norm();
    Eigen::Vector3d u1, v1 = a1 / a1mag + xhat;
    double v1mag = v1.norm();
    if(v1mag > 1e-6){
      u1 = v1 / v1mag;
    }else{
      u1 = yhat;
    }
    Eigen::Quaterniond q1(cos(M_PI/2.0), u1[0]*sin(M_PI/2.0), u1[1]*sin(M_PI/2.0), u1[2]*sin(M_PI/2.0));

    // Find quaternion to rotate second lattice vector to xy plane after applying above rotation
    Eigen::Vector3d a2prime = q1 * a2;
    double angle = -1.0 * atan2(a2prime[2], a2prime[1]);
    Eigen::Quaterniond q2(cos(angle/2.0), xhat[0]*sin(angle/2.0), xhat[1]*sin(angle/2.0), xhat[2]*sin(angle/2.0));

    // Calculate the new box
    double Lx = a1.norm();
    double a2x = a1.dot(a2) / Lx;
    double Ly = sqrt(a2.dot(a2) - a2x*a2x);
    double xy = a2x / Ly;
    Eigen::Vector3d v0xv1 = a1.cross(a2);
    double v0xv1mag = v0xv1.norm();
    double Lz = a3.dot(v0xv1) / v0xv1mag;
    double a3x = a1.dot(a3) / Lx;
    double xz = a3x / Lz;
    double yz = (a2.dot(a3) - a2x*a3x) / (Ly*Lz);

    box_ << Lx, xy*Ly, xz*Lz,
            0,  Ly,    yz*Lz,
            0,  0,     Lz;

    // This quaternion should rotate all positions into the new box
    Eigen::Quaterniond q = q2 * q1;
    for(size_t i = 0; i < positions_.size(); i++){
      positions_[i] = q * positions_[i];
    }
  }

  // Calculates the associated Legendre polynomial P_l^m(x) for positive m
  double MinkowskiStructureCalculator::LegendrePlm_m_gtr_0(int l, int m, double x)const{
    // Plm(x) = (-1^m / (2^l l!)) * (1-x^2)^(m/2) * (d/dx)^(l+m) (x^2-1)^l
    // A clever way to calculate Legendre polynomials based on Michiel Hermes' bond order code.
    if (m < 0 || m > l || fabs(x) > 1.0){
      std::cerr << "Bad arguments in MSM:LegendrePlm"<<l<<' '<<m<<' '<<x<<'\n';
      exit(42);
    }
    // Calculate the case where l==m: Pmm(x)
    double pmm = 1.0; // P00(x) = 1.0
    double sumx2 = sqrt(1-x*x);
    double fact = 1.0;
    for(int i = 1; i <= m; i++){
      pmm *= -fact * sumx2;
      fact += 2.0;
    }
    // (TeX notation) Go from P_m^m to P_l^m using two relations:
    // P_{l+1}^l(x) = x(2l+1)P_l^l(x)
    // (l-m)P_l^m(x) = x(2l+1)P_{l-1}^m(x) - (l+m-1)P_{l-2}^m(x)
    double pll = 0.0; // Only has an initial value to silence a compiler warning
    if(l != m){
      double pmmp1 = x*(2*m+1)*pmm; // P_{m+1}^m
      if(l == (m+1)){ return pmmp1; }
      else{
        for(int ll = m+2; ll <= l; ll++){
          pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm) / (ll-m);
          pmm = pmmp1;
          pmmp1 = pll;
        }
        return pll;
      }
    }
    return pmm;
  }
  // Unit test for the LegendrePlm_m_gtr_0 function
  void MinkowskiStructureCalculator::verify_LegendreP(void)const{
    double tolerance = 1e-14;
    // Reference values were obtained from Mathematica
    double ref_1 = 1.0;
    double ref_2 = -0.9927738916792684;
    double ref_3 = -41809.40924365284;
    double ref_4 = -38506.175717768485;
    double ref_5 = 1.0006646698482125e252;
    double test_1 = LegendrePlm_m_gtr_0(0, 0, 0.11);
    double test_2 = LegendrePlm_m_gtr_0(1, 1, 0.12);
    double test_3 = LegendrePlm_m_gtr_0(8, 6, 0.15);
    double test_4 = LegendrePlm_m_gtr_0(8, 6, -0.16);
    double test_5 = LegendrePlm_m_gtr_0(128, 128, -0.17);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_5, ref_5, tolerance, "unit test verify_LegendreP");
  }

  // Calculates ln((x-1)!) for integers x > 0 and tabulates them for reuse.
  double MinkowskiStructureCalculator::gammln_i(int x)const{
    // We tabulate the outcomes for efficiency
    static std::vector<double> table;
    if(x<1024){ // Don't build the table to unlimited size
      if( x+1 > int(table.size()) ){ // data for x hasn't been generated yet
        for(int xi = table.size(); xi <= x; xi++){ // so generate intermediate values up to x
          table.push_back( lgamma(double(xi)) );
        }
      }
      return table[x];
    } else return lgamma(double(x));
  }

  // Calculate the Wigner 3j-symbols using the Racah formula (up to l~32)
  double MinkowskiStructureCalculator::specialized_wigner3j(int l, int m1, int m2)const{
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    // Note: this is the 3j symbol only for ( l   l   l      )
    // Not the general one.                 ( m1  m2  -m1-m2 )
    // Tabulate calculated values for efficiency
    static std::map<std::tuple<int,int,int>, double> w3jMap;
    // Check if the Wigner 3j-symbol for this (l,m1,m2) is known already
    std::tuple<int,int,int> lm1m2_key = std::make_tuple(l,m1,m2);
    if(w3jMap.find(lm1m2_key) == w3jMap.end()){
      // If symbol is not known, calculate it and store it in the map.
      // But first, check if current algorithm can even handle requested l.
      if(l > 32){
        std::cerr << "Fatal error: cannot compute wigner3j symbols for l>32.";
        std::cerr << "Due to limited precision method, these would become inf/NaN.\n";
        exit(42);
      }
      // Okay, now calculate the wigner 3j-symbol.
      int m3 = -m1-m2;
      double wigner_factor = exp(3*gammln_i(l+1) - gammln_i(3*l+2)); // (l!)^3 / (3l+1)!
      // (-1)^m3 * sqrt(prefactor * (l-m1)!*(l+m1)!*(l-m2)!*(l+m2)!*(l-m3)!*(l-m3)! )
      double wigner = ((std::abs(m3) % 2 == 0)? 1: -1) * sqrt(wigner_factor * exp(
        gammln_i(l - m1 + 1) + gammln_i(l + m1 + 1) +
        gammln_i(l - m2 + 1) + gammln_i(l + m2 + 1) +
        gammln_i(l - m3 + 1) + gammln_i(l + m3 + 1)));
      double txsum = 0.0;
      for(int t = std::max(0, std::max(-m1,m2)); t <= std::min(l, std::min(l-m1,l+m2)); t++){
        // (-1)^t * t! * (l-t)! * (l-m1-t)! * (l+m2-t)! * (t+m1)! * (t-m2)!
        txsum += ((t % 2 == 0)? 1: -1) / exp(
          gammln_i(t + 1)          + gammln_i(l - t + 1)      +
          gammln_i(l - m1 - t + 1) + gammln_i(l + m2 - t + 1) +
          gammln_i(m1 + t + 1)     + gammln_i(-m2 + t + 1)   );
      }
      wigner *= txsum;
      if(w3jMap.size() < 1e6){ // Don't grow the map infinitely
        w3jMap[lm1m2_key] = wigner;
      } else return wigner;
    }
    return w3jMap[lm1m2_key];
  }
  // Unit test for the wigner3j function
  void MinkowskiStructureCalculator::verify_wigner3j(void)const{
    double tolerance = 1e-12;
    // Reference values were obtained from Mathematica
    double ref_1 = 1.0;                        // l=0,  m1=0,   m2=0
    double ref_2 = 1.0 / sqrt(6.0);            // l=1,  m1=0,   m2=1
    double ref_3 = -sqrt(3.0/35.0);            // l=2,  m1=2,   m2=-1
    double ref_4 = -3.0 / (2.0 * sqrt(143.0)); // l=5,  m1=-4,  m2=-1
    double ref_5 = 0.016988287171668218048;    // l=32, m1=-31, m2=17
    double test_1 = specialized_wigner3j(0, 0, 0);
    double test_2 = specialized_wigner3j(1, 0, 1);
    double test_3 = specialized_wigner3j(2, 2, -1);
    double test_4 = specialized_wigner3j(5, -4, -1);
    double test_5 = specialized_wigner3j(32, -31, 17);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_5, ref_5, tolerance, "unit test verify_wigner3j");
  }

  // Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
  double MinkowskiStructureCalculator::spherical_harmonic_factor(unsigned int l, int m)const{
    // We tabulate the constant prefactors for efficiency
    static std::vector<std::vector<double>> factor_table;
    if( l+1 > factor_table.size() ){ // data for l hasn't been generated yet
      for(size_t l_idx = 0; l_idx <= l; l_idx++){ // so generate prefactors up to l
        if( l_idx+1 > factor_table.size() ){ // only generate data that hasn't been generated already
          // std::cout << "l = " << l_idx << '\n';
          factor_table.resize(l_idx+1);
          factor_table[l_idx].resize(2*l_idx + 1);
          // std::cout << factor_table.size() << '\n';
          // std::cout << factor_table[l_idx].size() << '\n';
          for(size_t idx = 0; idx < factor_table[l_idx].size(); idx++){
            int m_idx = idx - l_idx; // first element of table is m = -l
            factor_table[l_idx][idx] = exp( gammln_i(l_idx - m_idx + 1) - gammln_i(l_idx + m_idx + 1) );
            // std::cout << " m = " << m_idx << ": " << factor_table[l_idx][idx] << '\n';
          }
        }
      }
    }
    return factor_table[l][m+l]; // offset by l because we can't have negative array indices
  }

  // Calculates a spherical harmonic Y_l^m(z,phi)
  std::complex<double> MinkowskiStructureCalculator::spherical_harmonic(int l, int m, double theta, double phi)const{
    // For the Minkowski metrics, the term "lfactor" cancels out in the final expression
    // for q_l. We retain it here to keep the expression for the spherical harmonics
    // complete and self-contained.
    std::complex<double> Ylm;
    double lfactor = (2*l + 1) / (4 * M_PI);
    // std::cerr << "Ylm " << l << ' ' << m << '\n';
    double lmfactorial = spherical_harmonic_factor(l, m);
    // The associated Legendre polynomials for positive and negative m are defined differently
    if(m >= 0){
      double Legendre_polynomial = sqrt(lfactor*lmfactorial) * LegendrePlm_m_gtr_0(l, m, cos(theta));
      Ylm.real( Legendre_polynomial * cos(m * phi) );
      Ylm.imag( Legendre_polynomial * sin(m * phi) );
    }else{
      double Legendre_polynomial = pow(-1.0,-m) * sqrt(lfactor/lmfactorial) * LegendrePlm_m_gtr_0(l, -m, cos(theta));
      Ylm.real( Legendre_polynomial * cos(m * phi) );
      Ylm.imag( Legendre_polynomial * sin(m * phi) );
    }
    return Ylm;
  }
  // Unit check for the spherical_harmonic function
  void MinkowskiStructureCalculator::verify_spherical_harmonic(void)const{
    double tolerance = 1e-5;
    // Reference values were obtained from Mathematica
    std::complex<double> ref_1(0.282095, 0.0);             // l=0, m=0,  theta=0.123, phi=0.321
    std::complex<double> ref_2(0.71284, 0.0);              // l=3, m=0,  theta=0.123, phi=0.321
    std::complex<double> ref_3(0.0122278,0.00914218);      // l=3, m=2,  theta=0.123, phi=0.321
    std::complex<double> ref_4(0.0122278,-0.00914218);     // l=3, m=-2, theta=0.123, phi=0.321
    std::complex<double> ref_5(9.03322e-06,-0.000263998);  // l=8, m=5,  theta=0.123, phi=0.321
    std::complex<double> ref_6(-9.03322e-06,-0.000263998); // l=8, m=-5, theta=0.123, phi=0.321
    std::complex<double> test_1 = spherical_harmonic(0, 0, 0.123, 0.321);
    std::complex<double> test_2 = spherical_harmonic(3, 0, 0.123, 0.321);
    std::complex<double> test_3 = spherical_harmonic(3, 2, 0.123, 0.321);
    std::complex<double> test_4 = spherical_harmonic(3, -2, 0.123, 0.321);
    std::complex<double> test_5 = spherical_harmonic(8, 5, 0.123, 0.321);
    std::complex<double> test_6 = spherical_harmonic(8, -5, 0.123, 0.321);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_5, ref_5, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_6, ref_6, tolerance, "unit test verify_spherical_harmonic");
  }

  // Does a number of checks to ensure the Voro++ output is usable
  void MinkowskiStructureCalculator::verify_voro_results(void)const{
    size_t N = pData_.size();
    // Verify that the number of neighbours equals the number of facets
    for(size_t i = 0; i < N; i++){
      if(pData_[i].nb_face_areas.size() != pData_[i].nb_indices.size()){
        std::cerr << "Fatal error: Voro++ number of neighbours does not match number of facets!\n";
        exit(42);
      }
    }
    // Verify that Voro neighbours are symmetric (if i nbs j, j nbs i), as they should be
    for(int i = 0; i < int(N); i++){
      for(int nb_i : pData_[i].nb_indices){
        auto found = std::find(
          std::begin(pData_[nb_i].nb_indices),
          std::end(pData_[nb_i].nb_indices),
          i
        );
        if(found == pData_[i].nb_indices.end()){
          std::cerr << "Fatal error: Voro++ neighbours are not symmetric!\n";
          exit(42);
        }
      }
    }
    // Verify that the sum of face areas matches the total cell area
    for(size_t i = 0; i < N; i++){
      double total_face_area = 0;
      for(double area : pData_[i].nb_face_areas){
        total_face_area += area;
      }
      if(fabs(total_face_area - pData_[i].total_face_area) / pData_[i].total_face_area > 1e-5){
        std::cerr << "Fatal error: Voro++ face area does not match." << '\n';
        std::cerr << "Total area of faces is " << total_face_area;
        std::cerr << " while cell area is " << pData_[i].total_face_area << '\n';
        exit(42);
      }
    }
    // Verify that Voro neighbours are unique. This can be false for very small systems.
    for(int i = 0; i < int(N); i++){
      for(size_t idx = 0; idx < pData_[i].nb_indices.size(); ++idx){
        int nb = pData_[i].nb_indices[idx];
        for(size_t idx_2 = 0; idx_2 < pData_[i].nb_indices.size(); ++idx_2){
          int nb_2 = pData_[i].nb_indices[idx_2];
          if(idx != idx_2 && nb == nb_2){
            std::cerr << "Fatal error: Voro++ neighbours are not unique!\n";
            std::cerr << "This can happen for very small systems, but bond order parameters will be wrong.\n";
            exit(42);
          }
        }
      }
    }
    // Verify that we're not looping around the periodic box
    for(size_t i = 0; i < N; i++){
      // Voro++ sets the particle as its own neighbour if its Voronoi cell
      // percolates the periodic volume. So we check for that.
      for(size_t j = 0; j < pData_[i].nb_indices.size(); j++){
        if(pData_[i].nb_indices[j] == int(i)){
          std::cerr << "Fatal error: Voro++ cell neighbours itself! Exiting.\n";
          exit(42);
        }
      }
    }
  }

  // Checks whether the current system is too small to yield valid Voro++ output
  bool MinkowskiStructureCalculator::is_system_too_small(void)const{
    // Verify that Voro neighbours are unique. This can be false for very small systems.
    for(int i = 0; i < int(pData_.size()); i++){
      for(size_t idx = 0; idx < pData_[i].nb_indices.size(); ++idx){
        int nb = pData_[i].nb_indices[idx];
        for(size_t idx_2 = 0; idx_2 < pData_[i].nb_indices.size(); ++idx_2){
          int nb_2 = pData_[i].nb_indices[idx_2];
          if(idx != idx_2 && nb == nb_2){
            std::cerr << "Warning: Voro++ neighbours are not unique!\n";
            std::cerr << "If your input system is very small (e.g. a unit cell) ";
            std::cerr << "you can probably safely ignore this warning.\n";
            return true;
          }
        }
      }
    }
    return false;
  }

  // Copies the input configuration a few times if it is too small
  void MinkowskiStructureCalculator::enlarge_system(void){
    // We enlarge the system by copying in on both sides of the original in x, y and z
    Eigen::Matrix3d newBox = 3 * box_ ;
    // Create a new set of positions by copying positions offset by the lattice vectors
    Eigen::Vector3d a1 = box_.col(0), a2 = box_.col(1), a3 = box_.col(2);
    std::vector<Eigen::Vector3d> newPositions = positions_;
    for(size_t idx = 0; idx < positions_.size(); idx++){
      // Then make new positions by copying the unit cell
      for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
          for(int k = -1; k <= 1; k++){
            if(i==0 && j==0 && k==0){ continue; } // Skip position of particle itself
            newPositions.push_back(positions_[idx] + i*a1 + j*a2 + k*a3);
          }
        }
      }
    }
    // Offset new positions to put them back into the positive xyz octant
    Eigen::Vector3d offset = box_.col(0) + box_.col(1) + box_.col(2);
    for(Eigen::Vector3d pos : newPositions){ pos += offset; }
    // Make these the new box and positions
    box_ = newBox;
    positions_ = newPositions;
  }

  // Loads a configuration and calculates the neighbour information using Voro++
  void MinkowskiStructureCalculator::load_configuration(
    const std::vector<std::vector<double>>& positions,
    const std::vector<double>& a1,
    const std::vector<double>& a2,
    const std::vector<double>& a3
  ){
    // Store the lattice vectors a1, a2, a3 in the columns of a box matrix
    box_ << a1[0], a2[0], a3[0],
            a1[1], a2[1], a3[1],
            a1[2], a2[2], a3[2];
    // Copy the input positions to the internal storage
    positions_.resize( positions.size() );
    for(size_t i = 0; i < positions_.size(); i++){
      positions_[i] << positions[i][0], positions[i][1], positions[i][2];
    }
    // Voro++ needs a1 to be along +x and a2 along +xy, so rotate if that's not the case
    if(a1[1] != 0 || a1[2] != 0 || a2[2] != 0){
      rotate_box_to_uppertriangular();
    }
    // Create a Voro++ box that will handle the Voronoi construction
    con_ = new voro::container_periodic(box_(0),
                                        box_(3), box_(4),
                                        box_(6), box_(7), box_(8),
                                        1, 1, 1,
                                        8);
    // Pass positions to Voro++
    for(size_t i = 0; i < positions_.size(); i++){
      con_->put(i, positions_[i][0], positions_[i][1], positions_[i][2]);
    }
    // Obtain neighbours from and areas and normals of the facets of the Voronoi cells
    pData_.resize( positions_.size() );
    voro::c_loop_all_periodic cloop(*con_);
    voro::voronoicell_neighbor c;
    if( cloop.start() ) do if( con_->compute_cell(c,cloop) ){
      c.neighbors( pData_[cloop.pid()].nb_indices   );
      c.face_areas(pData_[cloop.pid()].nb_face_areas);
      pData_[cloop.pid()].total_face_area = c.surface_area();
    } while(cloop.inc());

    // If the system is too small, Voro++ will yield unusable output in the
    // form of duplicate neighbours. To avoid this, we enlarge the system by
    // copying it until it is large enough.
    size_t copies = 0;
    while( is_system_too_small() ){
      if(copies > 3){
        std::cerr << "Fatal error: could not enlarge system enough to be valid.\n";
        exit(42);
      }
      enlarge_system();
      pData_.resize( positions_.size() );
      copies++;
      delete con_;
      con_ = new voro::container_periodic(box_(0),
                                          box_(3), box_(4),
                                          box_(6), box_(7), box_(8),
                                          1, 1, 1,
                                          8);
      for(size_t i = 0; i < positions_.size(); i++){
        con_->put(i, positions_[i][0], positions_[i][1], positions_[i][2]);
      }
      // Obtain neighbours from and areas and normals of the facets of the Voronoi cells
      pData_.resize( positions_.size() );
      voro::c_loop_all_periodic cloop(*con_);
      voro::voronoicell_neighbor c;
      if( cloop.start() ) do if( con_->compute_cell(c,cloop) ){
        c.neighbors( pData_[cloop.pid()].nb_indices   );
        c.face_areas(pData_[cloop.pid()].nb_face_areas);
        pData_[cloop.pid()].total_face_area = c.surface_area();
      } while(cloop.inc());
    }

    // Verify that the Voro++ output is now sensible
    verify_voro_results();

    // Create the bond angle vectors
    double dx,dy, dz, theta, phi;
    for(size_t i = 0; i < positions_.size(); i++){
      // for(size_t j = 0; j < pData_[i].nb_indices.size(); j++){
      for(int nb_i : pData_[i].nb_indices){
        Eigen::Vector3d& pos1 = positions_[i];
        Eigen::Vector3d& pos2 = positions_[nb_i];
        Eigen::Vector3d dr = nearest_image(pos1, pos2, box_);
        dx = dr[0];
        dy = dr[1];
        dz = dr[2];
        // theta (polar) [0,pi], phi (azimuthal) [-pi,pi]
        if(fabs(dx) < 1e-6 && fabs(dy) < 1e-6){ theta = dz > 0.0 ? 0 : M_PI; } // handle acos(1)
        else{ theta = acos(dz / sqrt(dx*dx + dy*dy + dz*dz)); }
        phi = atan2(dy,dx);
        // printf("%f %f %f %f %f\n",dx,dy,dz,theta,phi);
        pData_[i].thetas.push_back(theta);
        pData_[i].phis.push_back(phi);
      }
    }

    // Clear any stored qlm data from previously loaded configurations
    all_qlms_.clear();
    all_qlm_avs_.clear();
  }

  // Compute structure metrics for single particles and a single l.
  // These do not store their results, and should only be used if you
  // only need their values for a few particles in a larger system.
  std::complex<double> MinkowskiStructureCalculator::qlm(size_t p_idx, unsigned int l, int m)const{
    const particleData& p = pData_[p_idx];
    std::complex<double> qlm(0.0, 0.0);
    // Loop over the neighbours / facets
    for(size_t j = 0; j < p.nb_indices.size(); j++){
      // Calculate the weight factor from the area contribution
      double area_weight = p.nb_face_areas[j] / p.total_face_area;
      // Calculate the spherical harmonic
      std::complex<double> Ylm = spherical_harmonic(l, m, p.thetas[j], p.phis[j]);
      qlm += area_weight * Ylm;
    }
    return qlm;
  }
  std::complex<double> MinkowskiStructureCalculator::qlm_av(size_t p_idx, unsigned int l, int m)const{
    const particleData& p = pData_[p_idx];
    // Average over the neighbours, including itself
    std::complex<double> qlm_av = qlm(p_idx, l, m);
    for(size_t j = 0; j < p.nb_indices.size(); j++){
      qlm_av += qlm(j, l, m);
    }
    // Normalize average by number of neighbours
    return qlm_av / double(p.nb_indices.size()+1);
  }
  double MinkowskiStructureCalculator::ql(size_t p_idx, unsigned int l)const{
    double sum_m = 0.0;
    double factor = 4 * M_PI / (2*l + 1);
    for(int m = -int(l); m <= int(l); m++){
      sum_m += norm( qlm(p_idx, l, m) );
    }
    return sqrt(factor * sum_m);
  }
  double MinkowskiStructureCalculator::wl(size_t p_idx, unsigned int l)const{
    // Need the qlms first
    std::complex<double> qlms[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms[l+m] = qlm(p_idx, l, m);
    }

    // Now we can calculate the wl
    double wl = 0;
    // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    for(int m1 = -int(l); m1 <= int(l); m1++){
      for(int m2 = -int(l); m2 <= int(l); m2++){
        int m3 = -m1-m2;
        if(m3 < -int(l) || m3 > int(l)){continue;} // enforce -l <= m3 <= l
        wl += specialized_wigner3j(l,m1,m2) * (qlms[l+m1] * qlms[l+m2] * qlms[l+m3]).real();
      }
    }
    // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
    double qlms_norm = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      qlms_norm += norm(qlms[l+m]);
    }
    qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
    // Prevent errors if q's are almost 0 (e.g. for q=1, which should always be 0)
    if(wl < 1e-6 && qlms_norm < 1e-6){ return 0.0; }
    return wl / qlms_norm;
  }
  double MinkowskiStructureCalculator::ql_av(size_t p_idx, unsigned int l)const{
    double sum_m = 0.0;
    double factor = 4 * M_PI / (2*l + 1);
    for(int m = -int(l); m <= int(l); m++){
      sum_m += norm( qlm_av(p_idx, l, m) );
    }
    return sqrt(factor * sum_m);
  }
  double MinkowskiStructureCalculator::wl_av(size_t p_idx, unsigned int l)const{
    // Need the qlms first
    std::complex<double> qlms[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms[l+m] = qlm_av(p_idx, l, m);
    }

    // Now we can calculate the wl
    double wl = 0;
    // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    for(int m1 = -int(l); m1 <= int(l); m1++){
      for(int m2 = -int(l); m2 <= int(l); m2++){
        int m3 = -m1-m2;
        if(m3 < -int(l) || m3 > int(l)){continue;} // enforce -l <= m3 <= l
        wl += specialized_wigner3j(l,m1,m2) * (qlms[l+m1] * qlms[l+m2] * qlms[l+m3]).real();
      }
    }
    // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
    double qlms_norm = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      qlms_norm += norm(qlms[l+m]);
    }
    qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
    // Prevent errors if q's are almost 0 (e.g. for q=1, which should always be 0)
    if(wl < 1e-6 && qlms_norm < 1e-6){ return 0.0; }
    return wl / qlms_norm;
  }


  // These functions use all_qlms_ and all_qlm_avs_ to make sure we only compute
  // the qlm / avg qlm for each particle once. This is a significant optimization
  // when computing the averaged q/w's, or even when computing both q's and w's.
  const std::vector<std::complex<double>>& MinkowskiStructureCalculator::qlm_all(unsigned int l, int m){
    // Check if the qlms for this (l,m) pair are known already
    std::pair<unsigned int, int> lm_key = std::make_pair(l,m);
    if( all_qlms_[lm_key].empty() ){
      // If they are not known, calculate and store them
      std::vector<std::complex<double>> qlms( pData_.size(), std::complex<double>(0.0, 0.0) );
      for(size_t i = 0; i < pData_.size(); i++){
        qlms[i] = qlm(i, l, m);
      }
      all_qlms_[lm_key] = qlms;
    }
    return all_qlms_[lm_key];
  }
  std::vector<double> MinkowskiStructureCalculator::ql_all(unsigned int l){
    std::vector<double> qls(pData_.size(), 0.0);
    double factor = 4 * M_PI / (2*l + 1);
    for(size_t i = 0; i < pData_.size(); i++){
      double sum_m = 0.0;
      for(int m = -int(l); m <= int(l); m++){
        std::complex<double> qlm = qlm_all(l,m)[i];
        sum_m += norm(qlm);
      }
      qls[i] = sqrt(factor * sum_m);
    }
    return qls;
  }
  std::vector<double> MinkowskiStructureCalculator::wl_all(unsigned int l){
    // Get the required qlms in an (l, p_idx) structured array
    std::vector<std::complex<double>> qlms[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms[l+m] = qlm_all(l, m);
    }

    // Then calculate the wl's
    std::vector<double> wls(pData_.size());
    for(size_t i = 0; i < pData_.size(); i++){
      double wl = 0.0;
      // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
      // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
      for(int m1 = -int(l); m1 <= int(l); m1++){
        for(int m2 = -int(l); m2 <= int(l); m2++){
          int m3 = -m1-m2;
          if(m3 < -int(l) || m3 > int(l)){continue;} // enforce -l <= m3 <= l
          wl += specialized_wigner3j(l,m1,m2) * (qlms[l+m1][i] * qlms[l+m2][i] * qlms[l+m3][i]).real();
        }
      }
      // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
      double qlms_norm = 0.0;
      for(int m = -int(l); m <= int(l); m++){
        qlms_norm += norm(qlms[l+m][i]);
      }
      qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
      // Prevent errors if q's are almost 0 (e.g. for q=1, which should always be 0)
      if(wl < 1e-6 && qlms_norm < 1e-6){ wls[i] = 0.0; }
      else{ wls[i] = wl / qlms_norm; }
    }
    return wls;
  }
  const std::vector<std::complex<double>>& MinkowskiStructureCalculator::qlm_av_all(unsigned int l, int m){
    // Check if the averaged qlms for this (l,m) pair are known already
    std::pair<unsigned int, int> lm_key = std::make_pair(l,m);
    if( all_qlm_avs_[lm_key].empty() ){
      // If not, calculate them.
      // First, make sure all qlms for this (l,m) are available
      const std::vector<std::complex<double>>& qlms = qlm_all(l,m);
      // The compute qlm_av by averaging over neighbours plus itself
      std::vector<std::complex<double>> qlm_avs( pData_.size(), std::complex<double>(0.0, 0.0) );
      for(size_t i = 0; i < pData_.size(); i++){
        std::complex<double> qlm_av = qlms[i];
        for(size_t nb_i : pData_[i].nb_indices){
          qlm_av += qlms[nb_i];
        }
        qlm_avs[i] = qlm_av / double(pData_[i].nb_indices.size()+1);
      }
      all_qlm_avs_[lm_key] = qlm_avs;
    }
    return all_qlm_avs_[lm_key];
  }
  std::vector<double> MinkowskiStructureCalculator::ql_av_all(unsigned int l){
    std::vector<double> ql_avs(pData_.size(), 0.0);
    double factor = 4 * M_PI / (2*l + 1);
    for(size_t i = 0; i < pData_.size(); i++){
      double sum_m = 0.0;
      for(int m = -int(l); m <= int(l); m++){
        std::complex<double> qlm_av = qlm_av_all(l,m)[i];
        sum_m += norm(qlm_av);
      }
      ql_avs[i] = sqrt(factor * sum_m);
    }
    return ql_avs;
  }
  std::vector<double> MinkowskiStructureCalculator::wl_av_all(unsigned int l){
    // Get the required qlms in an (l, p_idx) structured array
    std::vector<std::complex<double>> qlm_avs[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlm_avs[l+m] = qlm_av_all(l, m);
    }

    // Then calculate the wl's
    std::vector<double> wl_avs(pData_.size());
    for(size_t i = 0; i < pData_.size(); i++){
      double wl_av = 0.0;
      // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
      // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
      for(int m1 = -int(l); m1 <= int(l); m1++){
        for(int m2 = -int(l); m2 <= int(l); m2++){
          int m3 = -m1-m2;
          if(m3 < -int(l) || m3 > int(l)){continue;} // enforce -l <= m3 <= l
          wl_av += specialized_wigner3j(l,m1,m2) * (qlm_avs[l+m1][i] * qlm_avs[l+m2][i] * qlm_avs[l+m3][i]).real();
        }
      }
      // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
      double qlm_avs_norm = 0.0;
      for(int m = -int(l); m <= int(l); m++){
        qlm_avs_norm += norm(qlm_avs[l+m][i]);
      }
      qlm_avs_norm = sqrt(qlm_avs_norm * qlm_avs_norm * qlm_avs_norm);
      // Prevent errors if q's are almost 0 (e.g. for q=1, which should always be 0)
      if(wl_av < 1e-6 && qlm_avs_norm < 1e-6){ wl_avs[i] = 0.0; }
      else{ wl_avs[i] = wl_av / qlm_avs_norm; }
    }
    return wl_avs;
  }

  // Another kind of bond order parameter from DOI: 10.1039/C8CP05248D
  std::vector<double> MinkowskiStructureCalculator::ql_dot_all(size_t i, unsigned int l){
    // Based on Eslami et al. (DOI: 10.1039/C8CP05248D)
    std::vector<double> qls(pData_.size(), 0.0);
    // Get the required qlms in an (l, p_idx) structured array
    std::complex<double> qlms_i[2 * l + 1], qlms_j[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms_i[l+m] = qlm_all(l, m)[i];
    }
    // Calculate the norms
    double norm_i = 0.0, norm_j = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      norm_i += norm(qlms_i[l+m]);
    }
    norm_i = sqrt(norm_i);
    // Calculate the dot product
    std::complex<double> dot = 0.0;
    // Sum over neighbours
    for(size_t nb_i : pData_[i].nb_indices){
      // Sum over m
      for(int m = -int(l); m <= int(l); m++){
        qlms_j[l+m] = qlm_all(l, m)[j];
        // Prevent precision errors if qlm's are almost 0 (e.g. for q=1, which should always be 0)
        if( !(abs(qlms_i[l+m]) < 1e-6 && abs(norm_i) < 1e-6) ){
          qls[l][i] += (qlms_i[l+m] / norm_i) * conj(qlms_j[l+m] / norm_j);
        }
        norm_j += norm(qlms_j[l+m]);
      }
      norm_j = sqrt(norm_j);
    }
  }

  // Calculates the dot product of the qlms to define measure of crystallinity
  double MinkowskiStructureCalculator::bond_crystallinity(size_t i, size_t j, unsigned int l){
    // Based on Rein ten Wolde et al. (DOI: 10.1063/1.471721)
    // Get the required qlms in an (l, p_idx) structured array
    std::complex<double> qlms_i[2 * l + 1], qlms_j[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms_i[l+m] = qlm_all(l, m)[i];
      qlms_j[l+m] = qlm_all(l, m)[j];
    }
    // Calculate the norms
    double norm_i = 0.0, norm_j = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      norm_i += norm(qlms_i[l+m]);
      norm_j += norm(qlms_j[l+m]);
    }
    norm_i = sqrt(norm_i);
    norm_j = sqrt(norm_j);
    // Calculate the dot product
    std::complex<double> dot = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      // Prevent precision errors if qlm's are almost 0 (e.g. for q=1, which should always be 0)
      if( !(abs(qlms_i[l+m]) < 1e-6 && abs(norm_i) < 1e-6) ){
        dot += (qlms_i[l+m] / norm_i) * conj(qlms_j[l+m] / norm_j);
      }
    }
    // Result should be real
    if(dot.imag() > 1e-6){
      std::cerr << "Error: dot product of complex qlms has nonzero imaginary part.\n";
    }
    return dot.real();
  }

  // As above, but averages over all input l
  double MinkowskiStructureCalculator::bond_crystallinity_lav(size_t i, size_t j, std::vector<unsigned int> all_l){
    // Extension of Rein ten Wolde et al. (DOI: 10.1063/1.471721)
    double dot_av = 0.0;
    for(unsigned int l : all_l){
      dot_av += bond_crystallinity(i, j, l);
    }
    return dot_av / double(all_l.size());
  }


}
