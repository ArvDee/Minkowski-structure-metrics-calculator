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
    const float& value,
    const float& ref_value,
    const float tolerance,
    const std::string location
  ){
    float diff = value - ref_value;
    float relative_diff = diff / ref_value;
    if(relative_diff > tolerance){
      std::cerr << "Fatal error in Minkowski metric, location '"<<location<<"'\n";
      std::cerr << "Value, reference value, relative difference and tolerance are:\n";
      std::cerr << value <<' '<< ref_value <<' '<< relative_diff <<' '<< tolerance <<'\n';
      exit(42);
    }
  }
  // Overload function to also handle the complex values from the spherical harmonics
  void verifySimilarity(
    const std::complex<float>& value,
    const std::complex<float>& ref_value,
    const float tolerance,
    const std::string location
  ){
    std::complex<float> diff = value - ref_value;
    float relative_diff = abs(diff) / abs(ref_value);
    if(relative_diff > tolerance){
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
    delete con_;
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

	// Calculates the associated Legendre polynomial P_l^m(x) for positive m
	float MinkowskiStructureCalculator::LegendrePlm_m_gtr_0(int l, int m, float x)const{
		// Code copied from Michiel Hermes' bond order code. Thanks Michiel!
	  float fact,pll=0.0,pmm,pmmp1,somx2;
	  int i,ll;
	  if (m < 0 || m > l || fabs(x) > 1.0)
	    printf("Bad arguments in MSM:LegendrePlm %i %i %f\n",l,m,fabs(x));
	  pmm=1.0;
	  if (m > 0) {
	    somx2=sqrt((1.0-x)*(1.0+x));
	    fact=1.0;
	    for (i=1;i<=m;i++) {
	      pmm *= -fact*somx2;
	      fact += 2.0;
	    }
	  }
	  if (l == m)
	    return pmm;
	  else {
	    pmmp1=x*(2*m+1)*pmm;
	    if (l == (m+1))
	      return pmmp1;
	    else {
	      for (ll=m+2;ll<=l;ll++) {
	        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	        pmm=pmmp1;
	        pmmp1=pll;
	      }
	      return pll;
	    }
	  }
	}
  // Unit test for the LegendrePlm_m_gtr_0 function
  void MinkowskiStructureCalculator::verify_LegendreP(void)const{
    float tolerance = 1e-6;
    // Reference values were obtained from Mathematica
    float ref_1 = 1.0;
    float ref_2 = -9.92774e-1;
    float ref_3 = -4.18094e4;
    float ref_4 = -3.85062e4;
    float test_1 = LegendrePlm_m_gtr_0(0, 0, 0.11);
    float test_2 = LegendrePlm_m_gtr_0(1, 1, 0.12);
    float test_3 = LegendrePlm_m_gtr_0(8, 6, 0.15);
    float test_4 = LegendrePlm_m_gtr_0(8, 6, -0.16);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_LegendreP");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_LegendreP");
  }

  // Calculates ln((x-1)!) for integers x > 0 and tabulates them for reuse.
  float MinkowskiStructureCalculator::gammln_i(int x)const{
    // We tabulate the outcomes for efficiency
    static std::vector<float> table;
    if( x+1 > int(table.size()) ){ // data for x hasn't been generated yet
      for(int xi = table.size(); xi <= x; xi++){ // so generate intermediate values up to x
        table.push_back( lgamma(float(xi)) );
      }
    }
    return table[x];
  }

  // Calculate the Wigner 3j-symbols using the Racah formula
  float MinkowskiStructureCalculator::wigner3j(int l, int m1, int m2, int m3)const{
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    // Note: this is the 3j symbol only for ( l  l  l  )
    // Not the general one.                 ( m1 m2 m3 )
    float wigner_factor = exp(3*gammln_i(l+1) - gammln_i(3*l+2)); // (l!)^3 / (3l+1)!
    // (-1)^m3 * sqrt(prefactor * (l-m1)!*(l+m1)!*(l-m2)!*(l+m2)!*(l-m3)!*(l-m3)! )
    float wigner = ((std::abs(m3) % 2 == 0)? 1: -1) * sqrt(wigner_factor * exp(
      gammln_i(l - m1 + 1) + gammln_i(l + m1 + 1) +
      gammln_i(l - m2 + 1) + gammln_i(l + m2 + 1) +
      gammln_i(l - m3 + 1) + gammln_i(l + m3 + 1)));
    float txsum = 0.0;
    for(int t = std::max(0, std::max(-m1,m2)); t <= std::min(l, std::min(l-m1,l+m2)); t++){
      // (-1)^t * t! * (l-t)! * (l-m1-t)! * (l+m2-t)! * (t+m1)! * (t-m2)!
      txsum += ((t % 2 == 0)? 1: -1) / exp(
        gammln_i(t + 1)          + gammln_i(l - t + 1)      +
        gammln_i(l - m1 - t + 1) + gammln_i(l + m2 - t + 1) +
        gammln_i(m1 + t + 1)     + gammln_i(-m2 + t + 1)   );
    }
    wigner *= txsum;
    return wigner;
  }
  // Unit test for the wigner3j function
  void MinkowskiStructureCalculator::verify_wigner3j(void)const{
    float tolerance = 1e-16;
    // Reference values were obtained from Mathematica
    float ref_1 = 1.0;                        // l=0, m1=0,  m2=0,  m3=0
    float ref_2 = 1.0 / sqrt(6.0);            // l=1, m1=0,  m2=1,  m3=-1
    float ref_3 = -sqrt(3.0/35.0);            // l=2, m1=2,  m2=-1, m3=-1
    float ref_4 = -3.0 / (2.0 * sqrt(143.0)); // l=5, m1=-4, m2=-2, m3=5
    float test_1 = wigner3j(0, 0, 0, 0);
    float test_2 = wigner3j(1, 0, 1, -1);
    float test_3 = wigner3j(2, 2, -1, -1);
    float test_4 = wigner3j(5, -4, -2, 5);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_wigner3j");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_wigner3j");
  }

	// Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
	float MinkowskiStructureCalculator::spherical_harmonic_factor(unsigned int l, int m)const{
		// We tabulate the constant prefactors for efficiency
		static std::vector<std::vector<float>> factor_table;
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
	std::complex<float> MinkowskiStructureCalculator::spherical_harmonic(int l, int m, float theta, float phi)const{
		// For the Minkowski metrics, the term "lfactor" cancels out in the final expression
		// for q_l. We retain it here to keep the expression for the spherical harmonics
		// complete and self-contained.
		std::complex<float> Ylm;
		float lfactor = (2*l + 1) / (4 * M_PI);
    // std::cerr << "Ylm " << l << ' ' << m << '\n';
		float lmfactorial = spherical_harmonic_factor(l, m);
		// The associated Legendre polynomials for positive and negative m are defined differently
		if(m >= 0){
			float Legendre_polynomial = sqrt(lfactor*lmfactorial) * LegendrePlm_m_gtr_0(l, m, cos(theta));
			Ylm.real( Legendre_polynomial * cos(m * phi) );
			Ylm.imag( Legendre_polynomial * sin(m * phi) );
		}else{
			float Legendre_polynomial = pow(-1.0,-m) * sqrt(lfactor/lmfactorial) * LegendrePlm_m_gtr_0(l, -m, cos(theta));
			Ylm.real( Legendre_polynomial * cos(m * phi) );
			Ylm.imag( Legendre_polynomial * sin(m * phi) );
		}
		return Ylm;
	}
  // Unit check for the spherical_harmonic function
  void MinkowskiStructureCalculator::verify_spherical_harmonic(void)const{
    float tolerance = 1e-5;
    // Reference values were obtained from Mathematica
    std::complex<float> ref_1(0.282095, 0.0);             // l=0, m=0,  theta=0.123, phi=0.321
    std::complex<float> ref_2(0.71284, 0.0);              // l=3, m=0,  theta=0.123, phi=0.321
    std::complex<float> ref_3(0.0122278,0.00914218);      // l=3, m=2,  theta=0.123, phi=0.321
    std::complex<float> ref_4(0.0122278,-0.00914218);     // l=3, m=-2, theta=0.123, phi=0.321
    std::complex<float> ref_5(9.03322e-06,-0.000263998);  // l=8, m=5,  theta=0.123, phi=0.321
    std::complex<float> ref_6(-9.03322e-06,-0.000263998); // l=8, m=-5, theta=0.123, phi=0.321
    std::complex<float> test_1 = spherical_harmonic(0, 0, 0.123, 0.321);
    std::complex<float> test_2 = spherical_harmonic(3, 0, 0.123, 0.321);
    std::complex<float> test_3 = spherical_harmonic(3, 2, 0.123, 0.321);
    std::complex<float> test_4 = spherical_harmonic(3, -2, 0.123, 0.321);
    std::complex<float> test_5 = spherical_harmonic(8, 5, 0.123, 0.321);
    std::complex<float> test_6 = spherical_harmonic(8, -5, 0.123, 0.321);
    verifySimilarity(test_1, ref_1, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_2, ref_2, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_3, ref_3, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_4, ref_4, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_5, ref_5, tolerance, "unit test verify_spherical_harmonic");
    verifySimilarity(test_6, ref_6, tolerance, "unit test verify_spherical_harmonic");
  }

  // Does a number of checks to ensure the Voro++ output is usable
  void MinkowskiStructureCalculator::verify_voro_results(void)const{
    size_t N = nbs_.size();
    // Verify that the number of neighbours equals the number of facets
    for(size_t i = 0; i < N; i++){
      if(nbs_[i].face_areas.size() != nbs_[i].indices.size()){
        std::cerr << "Fatal error: Voro++ number of neighbours does not match number of facets!\n";
        exit(42);
      }
    }
    // Verify that Voro neighbours are symmetric (if i nbs j, j nbs i), as they should be
    for(int i = 0; i < int(N); i++){
      for(int nb_i : nbs_[i].indices){
        auto found = std::find(
          std::begin(nbs_[nb_i].indices),
          std::end(nbs_[nb_i].indices),
          i
        );
        if(found == nbs_[i].indices.end()){
          std::cerr << "Fatal error: Voro++ neighbours are not symmetric!\n";
          exit(42);
        }
      }
    }
    // Verify that the sum of face areas matches the total cell area
    for(size_t i = 0; i < N; i++){
      double total_face_area = 0;
      for(double area : nbs_[i].face_areas){
        total_face_area += area;
      }
      if(fabs(total_face_area - nbs_[i].cell_area) / nbs_[i].cell_area > 1e-5){
        std::cerr << "Fatal error: Voro++ face area does not match." << '\n';
        std::cerr << "Total area of faces is " << total_face_area;
        std::cerr << " while cell area is " << nbs_[i].cell_area << '\n';
        exit(42);
      }
    }
    // Verify that Voro neighbours are unique. This can be false for very small systems.
    for(int i = 0; i < int(N); i++){
      for(size_t idx = 0; idx < nbs_[i].indices.size(); ++idx){
        int nb = nbs_[i].indices[idx];
        for(size_t idx_2 = 0; idx_2 < nbs_[i].indices.size(); ++idx_2){
          int nb_2 = nbs_[i].indices[idx_2];
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
			for(size_t j = 0; j < nbs_[i].indices.size(); j++){
				if(nbs_[i].indices[j] == int(i)){
					std::cerr << "Fatal error: Voro++ cell neighbours itself! Exiting.\n";
					exit(42);
				}
			}
		}
  }

  // Checks whether the current system is too small to yield valid Voro++ output
  bool MinkowskiStructureCalculator::is_system_too_small(void)const{
    // Verify that Voro neighbours are unique. This can be false for very small systems.
    for(int i = 0; i < int(nbs_.size()); i++){
      for(size_t idx = 0; idx < nbs_[i].indices.size(); ++idx){
        int nb = nbs_[i].indices[idx];
        for(size_t idx_2 = 0; idx_2 < nbs_[i].indices.size(); ++idx_2){
          int nb_2 = nbs_[i].indices[idx_2];
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
  void MinkowskiStructureCalculator::msm_prepare(
    const std::vector<std::vector<float>>& positions,
    const std::vector<float>& box
  ){
    // Store the lattice vectors a1, a2, a3 in the columns of a box matrix
    box_ << box[0], box[3], box[6],
            box[1], box[4], box[7],
            box[2], box[5], box[8];
    // Voro++ needs a1 to be along +x, a2 along +xy and a3 in the +xyz octant.
    // TODO: write code to allow other boxes too by rotating them.
    if(box_(1) != 0 || box_(2) != 0 || box_(5) != 0){
      std::cerr << "MSM Error: Input box is not row-major upper triangular, which is required! Exiting.\n";
      exit(42);
    }
    // Create a Voro++ box that will handle the Voronoi construction
    con_ = new voro::container_periodic(box_(0),
																        box_(3), box_(4),
																        box_(6), box_(7), box_(8),
																        1, 1, 1,
																        8);
    // Copy the input positions to the internal storage and pass them to Voro++
    positions_.resize( positions.size() );
    for(size_t i = 0; i < positions_.size(); i++){
      positions_[i] << positions[i][0], positions[i][1], positions[i][2];
      con_->put(i, positions_[i][0], positions_[i][1], positions_[i][2]);
    }
    // Obtain neighbours from and areas and normals of the facets of the Voronoi cells
    nbs_.resize( positions_.size() );
    voro::c_loop_all_periodic cloop(*con_);
    voro::voronoicell_neighbor c;
    if( cloop.start() ) do if( con_->compute_cell(c,cloop) ){
      c.neighbors( nbs_[cloop.pid()].indices   );
      c.face_areas(nbs_[cloop.pid()].face_areas);
      nbs_[cloop.pid()].cell_area = c.surface_area();
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
      nbs_.resize( positions_.size() );
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
      nbs_.resize( positions_.size() );
      voro::c_loop_all_periodic cloop(*con_);
      voro::voronoicell_neighbor c;
      if( cloop.start() ) do if( con_->compute_cell(c,cloop) ){
        c.neighbors( nbs_[cloop.pid()].indices   );
        c.face_areas(nbs_[cloop.pid()].face_areas);
        nbs_[cloop.pid()].cell_area = c.surface_area();
      } while(cloop.inc());
    }

    // Verify that the Voro++ output is now sensible
    verify_voro_results();

    // Create the bond angle vectors
    float dx,dy, dz, theta, phi;
    for(size_t i = 0; i < positions_.size(); i++){
      // for(size_t j = 0; j < nbs_[i].indices.size(); j++){
      for(int nb_i : nbs_[i].indices){
        Eigen::Vector3d& pos1 = positions_[i];
        Eigen::Vector3d& pos2 = positions_[nb_i];
        Eigen::Vector3d dr = nearest_image(pos1, pos2, box_);
        dx = dr[0];
        dy = dr[1];
        dz = dr[2];
        // theta (polar) [0,pi], phi (azimuthal) [-pi,pi]
        if(dx == 0.0 && dy == 0.0){ theta = 0.0; } // handle acos(1)
        else{ theta = acos(dz / sqrt(dx*dx + dy*dy + dz*dz)); }
        phi = atan2(dy,dx);
        nbs_[i].thetas.push_back(theta);
        nbs_[i].phis.push_back(phi);
      }
    }
  }

  // Calculates q_l for particle p
  float MinkowskiStructureCalculator::ql(unsigned int p_idx, unsigned int l)const{
    float sum_m = 0.0;
		float factor = 4 * M_PI / (2*l + 1);
		for(int m = -int(l); m <= int(l); m++){
			std::complex<float> qlm = 0;
			// Loop over the neighbours / facets
			for(size_t j = 0; j < nbs_[p_idx].indices.size(); j++){
				// Calculate the weight factor from the area contribution
				float area_weight = nbs_[p_idx].face_areas[j] / nbs_[p_idx].cell_area;
				// Calculate the spherical harmonic
				std::complex<float> Ylm = spherical_harmonic(l, m, nbs_[p_idx].thetas[j], nbs_[p_idx].phis[j]);
				qlm += area_weight * Ylm;
        // if(l==1)printf("%d %d %f %f %f + %fi\n",l,m,nbs_[p_idx].thetas[j],nbs_[p_idx].phis[j],Ylm.real(),Ylm.imag());
			}
			sum_m += norm(qlm);
		}
		return sqrt(factor * sum_m);
	}

  // Calculates w_l for particle p
  float MinkowskiStructureCalculator::wl(unsigned int p_idx, unsigned int l)const{
    // Need the qlms first
    std::complex<float> qlms[2 * l + 1];
    for(int m = -int(l); m <= int(l); m++){
      qlms[l+m] = 0.0;
      // Loop over the neighbours / facets
      for(size_t j = 0; j < nbs_[p_idx].indices.size(); j++){
        // Calculate the weight factor from the area contribution
        float area_weight = nbs_[p_idx].face_areas[j] / nbs_[p_idx].cell_area;
        // Calculate the spherical harmonic
        std::complex<float> Ylm = spherical_harmonic(l, m, nbs_[p_idx].thetas[j], nbs_[p_idx].phis[j]);
        qlms[l+m] += area_weight * Ylm;
      }
    }

    // Now we can calculate the wl
    float wl = 0;
    // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    for(int m1 = -int(l); m1 <= int(l); m1++){
      for(int m2 = -int(l); m2 <= int(l); m2++){
        int m3 = -m1-m2;
        if(m3 < -int(l) || m3 > int(l)){continue;} // enforce -l <= m3 <= l
        wl += wigner3j(l,m1,m2,m3) * (qlms[l+m1] * qlms[l+m2] * qlms[l+m3]).real();
      }
    }
    // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
    float qlms_norm = 0.0;
    for(int m = -int(l); m <= int(l); m++){
      qlms_norm += norm(qlms[l+m]);
    }
    qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
    // Prevent errors if q's are almost 0 (e.g. for q=1, which should always be 0)
    if(wl < 1e-6 && qlms_norm < 1e-6){ return 0.0; }
    return wl / qlms_norm;
  }

  // Calculates q_l_av for particle p (q_l averaged over neighbours)
  // float MinkowskiStructureCalculator::ql_av(unsigned int p_idx, unsigned int l)const{
  //   // Need the qlms first
  //   std::complex<float> qlms[2 * l + 1];
  //   for(int m = -int(l); m <= int(l); m++){
  //     qlms[l+m] = 0.0;
  //     // Loop over the neighbours / facets
  //     for(size_t j = 0; j < nbs_[p_idx].indices.size(); j++){
  //       // Calculate the weight factor from the area contribution
  //       float area_weight = nbs_[p_idx].face_areas[j] / nbs_[p_idx].cell_area;
  //       // Calculate the spherical harmonic
  //       std::complex<float> Ylm = spherical_harmonic(l, m, nbs_[p_idx].thetas[j], nbs_[p_idx].phis[j]);
  //       qlms[l+m] += area_weight * Ylm;
  //     }
  //   }
  //   // Now calculate the averaged ql
  // }

  // Fills the q and w vectors with their values for input positions, starting from l=0
  void MinkowskiStructureCalculator::compute(
      const std::vector<std::vector<float>>& positions,
      const std::vector<float>& box,
      std::vector<std::vector<float>>& q,
      std::vector<std::vector<float>>& w
  ){
    msm_prepare(positions, box);
    // Calculate the q's
    for(size_t i = 0; i < q.size(); i++){
      for(size_t l = 0; l < q[i].size(); l++){
        q[i][l] = ql(i,l);
      }
    }
    // Calculate the w's
    for(size_t i = 0; i < w.size(); i++){
      for(size_t l = 0; l < w[i].size(); l++){
        w[i][l] = wl(i,l);
      }
    }
  }

}
