#include "MinkowskiStructureMetrics.h"

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


namespace MSM {

  Voro_Nbs nbs_;

	// Calculates the associated Legendre polynomial P_l^m(x) for positive m
	// RvD: Checked outputs, matches GNU scientific library with rel. error ~1E-8.
	float LegendrePlm_m_gtr_0(int l, int m, double x){
		// Code copied from Michiel Hermes' bond order code. Thanks Michiel!
	  double fact,pll=0.0,pmm,pmmp1,somx2;
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

	// Calculates ln((x-1)!) so we can calculate (l-m)!/(l+m)! as exp(gammln(l-m)+gammln(l+m))
	// RvD: Checked outputs, matches GNU scientific library with rel. error ~1E-16.
	double gammln(double xx){
		// Code copied from Michiel Hermes' bond order code. Thanks Michiel!
	  double x,y,tmp,ser;
	  static double cof[6]={76.18009172947146,    -86.50532032941677,
	                        24.01409824083091,    -1.231739572450155,
	                        0.1208650973866179e-2,-0.5395239384953e-5};
	  int j;
	  y=x=xx;
	  tmp=x+5.5;
	  tmp -= (x+0.5)*log(tmp);
	  ser=1.000000000190015;
	  for (j=0;j<=5;j++) ser += cof[j] / ++y;
	  return -tmp+log(2.5066282746310005*ser/x);
	}

	// Calculates prefactors (l-m)!/(l+m)! using log gamma functions for the factorials
	// RvD: outputs match Mathematica's output at least up to l=12.
	double spherical_harmonic_factor(unsigned int l, int m){
		// We tabulate the constant prefactors for efficiency
		static std::vector<std::vector<double>> factor_table;
		if( l+1 > factor_table.size() ){ // data for l hasn't been generated yet
			for(size_t l_idx = 0; l_idx <= l; l_idx++){ // so generate prefactors up to l
				if( l_idx+1 > factor_table.size() ){ // only generate data that hasn't been generated already
					// std::cout << "l = " << l << '\n';
					factor_table.resize(l_idx+1);
					factor_table[l_idx].resize(2*l + 1);
					// std::cout << factor_table.size() << '\n';
					// std::cout << factor_table[l_idx].size() << '\n';
					for(size_t idx = 0; idx < factor_table[l_idx].size(); idx++){
						int m_idx = idx - l_idx; // first element of table is m = -l
						factor_table[l_idx][idx] = exp( gammln(l_idx - m_idx + 1) - gammln(l_idx + m_idx + 1) );
						// std::cout << " m = " << m_idx << ": " << factor_table[l_idx][idx] << '\n';
					}
				}
			}
		}
		return factor_table[l][m+l]; // offset by l because we can't have negative array indices
	}

	// Calculates a spherical harmonic Y_l^m(z,phi)
	std::complex<double> spherical_harmonic(int l, int m, double theta, double phi){
		// For the Minkowski metrics, the term "lfactor" cancels out in the final expression
		// for q_l. We retain it here to keep the expression for the spherical harmonics
		// complete and self-contained.
		std::complex<double> Ylm;
		double lfactor = (2*l + 1) / (4 * M_PI);
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

	// C++ compatible format
	void msm(
		const std::vector<std::vector<double>>& positions,
		const std::vector<double>& box,
		std::vector<std::vector<double>>& q,
		std::vector<std::vector<double>>& w
	){
		size_t n_positions = positions.size();
    // A safety check
    if(q.size() != w.size() || q[0].size() != w[0].size()){
      std::cout << "Error: vectors for ql and wl are not the same size. Fix that first, please.\n";
      exit(42);
    }

		// Check whether box is column-major upper triangular (needed for voro++)
		if(box[3] != 0 || box[6] != 0 || box[7] != 0){
			std::cout << "MSM Error: Box is not row-major upper triangular, which is required! Exiting.\n";
			exit(42);
			// TODO: some automatic detection of other valid box inputs
			// e.g. row-major (don't forget to fix corresponding positions!)
		}

		// Define a box for the Voronoi tesselation
		voro::container_periodic con(box[0],
																 box[1], box[4],
																 box[2], box[5], box[8],
																 1, 1, 1, // number of internal subdivisions of box in voro++
																 8);
		// Add positions into the Voro++ box
		for(size_t i = 0; i < n_positions; i++){
			con.put(i, positions[i][0], positions[i][1], positions[i][2]);
		}
		// for (size_t i = 0; i < positions.size(); i++) {
		// 	std::cout << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << ' '<< '\n';
		// }

		// Obtain neighbours from and areas and normals of the facets of the Voronoi cells
		voro::c_loop_all_periodic cloop(con);
		voro::voronoicell_neighbor c;
		std::vector<std::vector<int> > neighbours(n_positions);
		std::vector<std::vector<double> > face_areas(n_positions);
		std::vector<double> cell_areas;
		if( cloop.start() ) do if( con.compute_cell(c,cloop) ){
			// cloop.pid() gives the index of the cell / particle currently being considered.
			c.neighbors(  neighbours[cloop.pid()] );
			c.face_areas( face_areas[cloop.pid()] );
			cell_areas.push_back( c.surface_area() );
		} while(cloop.inc());

		// for(size_t i = 0; i < neighbours[0].size(); i++){
		// 	std::cout << neighbours[0][i]
		// 	<<' '<< positions[ neighbours[0][i] ][0]
		//   <<' '<< positions[ neighbours[0][i] ][1]
		//   <<' '<< positions[ neighbours[0][i] ][2]
		// 	<< '\n';
		// }
		// std::cout << '\n';

		// Check that we're not looping around the periodic box
		for(size_t i = 0; i < q.size(); i++){
			// Voro++ sets the particle as its own neighbour if its Voronoi cell
			// percolates the periodic volume. So we check for that.
			for(size_t j = 0; j < neighbours[i].size(); j++){
				if(neighbours[i][j] == ((int) i)){
					std::cout << "Error: Voronoi cells neighbours itself! Exiting.\n";
					exit(42);
				}
			}
		}

		// We can now calculate the structure metrics using the Voronoi information
		double sum_m;
    // Loop over particles
		for(size_t i = 0; i < q.size(); i++){
			// Loop over all l that we want to know
      unsigned int n_l = q[i].size();
      std::complex<double> qlms[n_l][2 * n_l + 1];
			for(unsigned int l = 0; l < n_l; l++){
        // ql's
				sum_m = 0;
				double factor = 4 * M_PI / (2*l + 1);
				for(int m = -((int) l); m <= ((int) l); m++){
					qlms[l][l+m] = 0;
					// Loop over the neighbours / facets
					for(size_t j = 0; j < neighbours[i].size(); j++){
						// Calculate the weight factor from the area contribution
						double area_weight = face_areas[i][j] / cell_areas[i];
						// Get the angles of the bond
						double dx = positions[ neighbours[i][j] ][0] - positions[i][0];
						double dy = positions[ neighbours[i][j] ][1] - positions[i][1];
						double dz = positions[ neighbours[i][j] ][2] - positions[i][2];
						// theta (polar) [0,pi], phi (azimuthal) [-pi,pi]
						dz /= sqrt(dx*dx + dy*dy + dz*dz); // cos(theta) = dz/r
						double phi = atan2(dy,dx);
						// Calculate the spherical harmonic
						std::complex<double> Ylm = spherical_harmonic(l, m, acos(dz), phi); // TODO optimize angles
						qlms[l][l+m] += area_weight * Ylm;
					}
					sum_m += norm(qlms[l][l+m]);
				}
				q[i][l] = sqrt(factor * sum_m);
        // wl's
        double wigner_factor = exp(3*gammln(l+1) - gammln(3*l+2)); // (l!)^3 / (3l+1)!
        // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
        // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
        for(int m1 = -((int) l); m1 <= ((int) l); m1++){
          for(int m2 = -((int) l); m2 <= ((int) l); m2++){
            int m3 = -m1-m2;
            if(m3 < -((int) l) || m3 > ((int) l)){continue;} // enforce -l <= m3 <= l
            // (-1)^m3 * sqrt(prefactor * (l-m1)!*(l+m1)!*(l-m2)!*(l+m2)!*(l-m3)!*(l-m3)! )
            double wigner = ((std::abs(m3) % 2 == 0)? 1: -1) * sqrt(wigner_factor * exp(
              gammln(l - m1 + 1) + gammln(l + m1 + 1) +
              gammln(l - m2 + 1) + gammln(l + m2 + 1) +
              gammln(l - m3 + 1) + gammln(l + m3 + 1)));
            double txsum = 0.0;
            for(int t = std::max(0, std::max(-m1, m2)); t <= std::min(((int) l), std::min(-m1, m2) + ((int) l)); t++){
              // (-1)^t * t! * (l-t)! * (l-m1-t)! * (l+m2-t)! * (t+m1)! * (t-m2)!
              txsum += ((t % 2 == 0)? 1: -1) / exp(
                gammln(t + 1)          + gammln(l - t + 1)      +
                gammln(l - m1 - t + 1) + gammln(l + m2 - t + 1) +
                gammln(m1 + t + 1)     + gammln(-m2 + t + 1)   );
            }
            wigner *= txsum;
            // std::cout << "l m1 m2 w3j "<<l<<' '<<m1<<' '<<m2<<' '<<wigner<<'\n';
            w[i][l] += wigner * (qlms[l][l+m1] * qlms[l][l+m2] * qlms[l][l+m3]).real();
          }
        }
        // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
        double qlms_norm = 0.0;
				for(int m = -((int) l); m <= ((int) l); m++){
          qlms_norm += std::norm(qlms[l][m + l]);
        }
        qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
        w[i][l] /= qlms_norm;
      }
    }
	}

  // Loads a configuration and calculates the neighbour information using Voro++
  void msm_prepare(unsigned int n_positions, double **positions, double box[9]){
    // Erase existing neighbour data
    nbs_.indices.clear();
    nbs_.face_areas.clear();
    nbs_.cell_areas.clear();
    nbs_.thetas.clear();
    nbs_.phis.clear();
    // Check whether box is column-major upper triangular (needed for voro++)
    if(box[3] != 0 || box[6] != 0 || box[7] != 0){
      std::cout << "MSM Error: Box is not row-major upper triangular, which is required! Exiting.\n";
      exit(42);
    }
    // Create a Voro++ box
    voro::container_periodic con(box[0],
																 box[1], box[4],
																 box[2], box[5], box[8],
																 1, 1, 1, // number of internal subdivisions of box in voro++
																 8);
    // Fill it with the positions
    for(size_t i = 0; i < n_positions; i++){
			con.put(i, positions[i][0], positions[i][1], positions[i][2]);
		}
		// Obtain neighbours from and areas and normals of the facets of the Voronoi cells
		voro::c_loop_all_periodic cloop(con);
		voro::voronoicell_neighbor c;
		if( cloop.start() ) do if( con.compute_cell(c,cloop) ){
			c.neighbors(  nbs_.indices[   cloop.pid()] );
			c.face_areas( nbs_.face_areas[cloop.pid()] );
			nbs_.cell_areas.push_back( c.surface_area() );
		} while(cloop.inc());
		// Check that we're not looping around the periodic box
		for(size_t i = 0; i < n_positions; i++){
			// Voro++ sets the particle as its own neighbour if its Voronoi cell
			// percolates the periodic volume. So we check for that.
			for(size_t j = 0; j < nbs_.indices[i].size(); j++){
				if(nbs_.indices[i][j] == ((int) i)){
					std::cout << "Error: Voronoi cells neighbours itself! Exiting.\n";
					exit(42);
				}
			}
		}
    // Create the bond angle vectors
    for(size_t i = 0; i < n_positions; i++){
      for(size_t j = 0; j < nbs_.indices[i].size(); j++){
        double dx = positions[ nbs_.indices[i][j] ][0] - positions[i][0];
        double dy = positions[ nbs_.indices[i][j] ][1] - positions[i][1];
        double dz = positions[ nbs_.indices[i][j] ][2] - positions[i][2];
        // theta (polar) [0,pi], phi (azimuthal) [-pi,pi]
        double theta = acos(dz / sqrt(dx*dx + dy*dy + dz*dz)); // cos(theta) = dz/r
        double phi = atan2(dy,dx);
        nbs_.thetas[i].push_back(theta);
        nbs_.phis[i].push_back(phi);
      }
    }
  }

  // Calculates q_l for particle p
  double ql(unsigned int p_idx, unsigned int l){
    double sum_m = 0.0;
		double factor = 4 * M_PI / (2*l + 1);
		for(int m = -((int) l); m <= ((int) l); m++){
			std::complex<double> qlm = 0;
			// Loop over the neighbours / facets
			for(size_t j = 0; j < nbs_.indices[p_idx].size(); j++){
				// Calculate the weight factor from the area contribution
				double area_weight = nbs_.face_areas[p_idx][j] / nbs_.cell_areas[p_idx];
				// Calculate the spherical harmonic
				std::complex<double> Ylm = spherical_harmonic(l, m, nbs_.thetas[p_idx][j], nbs_.phis[p_idx][j]);
				qlm += area_weight * Ylm;
			}
			sum_m += norm(qlm);
		}
		return sqrt(factor * sum_m);
	}

  // Calculates w_l for particle p
  double wl(unsigned int p_idx, unsigned int l){
    // Need the qlms first
    std::complex<double> qlms[2 * l + 1];
    for(int m = -((int) l); m <= ((int) l); m++){
      qlms[l+m] = 0.0;
      // Loop over the neighbours / facets
      for(size_t j = 0; j < nbs_.indices[p_idx].size(); j++){
        // Calculate the weight factor from the area contribution
        double area_weight = nbs_.face_areas[p_idx][j] / nbs_.cell_areas[p_idx];
        // Calculate the spherical harmonic
        std::complex<double> Ylm = spherical_harmonic(l, m, nbs_.thetas[p_idx][j], nbs_.phis[p_idx][j]);
        qlms[l+m] += area_weight * Ylm;
      }
    }

    // Now we can calculate the wl
    double wl = 0;
    double wigner_factor = exp(3*gammln(l+1) - gammln(3*l+2)); // (l!)^3 / (3l+1)!
    // Calculate the wl, using the Racah formula for the Wigner 3j-symbols
    // (Quantum Mechanics Volume II, Albert Messiah, 1962, p.1058)
    for(int m1 = -((int) l); m1 <= ((int) l); m1++){
      for(int m2 = -((int) l); m2 <= ((int) l); m2++){
        int m3 = -m1-m2;
        if(m3 < -((int) l) || m3 > ((int) l)){continue;} // enforce -l <= m3 <= l
        // (-1)^m3 * sqrt(prefactor * (l-m1)!*(l+m1)!*(l-m2)!*(l+m2)!*(l-m3)!*(l-m3)! )
        double wigner = ((std::abs(m3) % 2 == 0)? 1: -1) * sqrt(wigner_factor * exp(
          gammln(l - m1 + 1) + gammln(l + m1 + 1) +
          gammln(l - m2 + 1) + gammln(l + m2 + 1) +
          gammln(l - m3 + 1) + gammln(l + m3 + 1)));
        double txsum = 0.0;
        for(int t = std::max(0, std::max(-m1, m2)); t <= std::min(((int) l), std::min(-m1, m2) + ((int) l)); t++){
          // (-1)^t * t! * (l-t)! * (l-m1-t)! * (l+m2-t)! * (t+m1)! * (t-m2)!
          txsum += ((t % 2 == 0)? 1: -1) / exp(
            gammln(t + 1)          + gammln(l - t + 1)      +
            gammln(l - m1 - t + 1) + gammln(l + m2 - t + 1) +
            gammln(m1 + t + 1)     + gammln(-m2 + t + 1)   );
        }
        wigner *= txsum;
        // std::cout << "l m1 m2 w3j "<<l<<' '<<m1<<' '<<m2<<' '<<wigner<<'\n';
        wl += wigner * (qlms[l+m1] * qlms[l+m2] * qlms[l+m3]).real();
      }
    }
    // Normalize wl by 1.0 / (|ql|^2)^(3/2) to map it into the range [0,1]
    double qlms_norm = 0.0;
    for(int m = -((int) l); m <= ((int) l); m++){
      qlms_norm += std::norm(qlms[l+m]);
    }
    qlms_norm = sqrt(qlms_norm * qlms_norm * qlms_norm);
    return wl / qlms_norm;
  }



    // input: set of coordinates in a box of which we want the Minkowski structure metrics
    // output: array of requested MSM

  	// c compatible format NOTE WIP
  	// void msmc(size_t n_positions, double **positions, double box[9], double (*q)[6]){
    //
  	// 	// Check whether box is column-major upper triangular (needed for voro++)
  	// 	if(box[3] != 0 || box[6] != 0 || box[7] != 0){
  	// 		std::cout << "MSM Error: Box is not row-major upper triangular, which is required! Exiting.\n";
  	// 		exit(42);
  	// 		// TODO: some automatic detection of other valid box inputs
  	// 		// e.g. column-major (don't forget to fix corresponding positions!)
  	// 	}
    //
  	// 	// Define a box for the Voronoi tesselation
  	// 	// TODO: add nx,ny,nz calculation
  	// 	voro::container_periodic con(box[0],
  	// 		 													 box[1], box[4],
  	// 															 box[2], box[5], box[8],
  	// 															 // nx,ny,nz, // number of internal subdivisions of box in voro++
  	// 															 1, 1, 1, // number of internal subdivisions of box in voro++
  	// 															 8);
  	// 	// Add positions into the Voro++ box
  	// 	for(size_t i = 0; i < n_positions; i++){
  	// 		con.put(i, positions[i][0], positions[i][1], positions[i][2]);
  	// 	}
    //
  	// 	// Obtain neighbours from and areas and normals of the facets of the Voronoi cells
  	// 	voro::c_loop_all_periodic cloop(con);
  	// 	voro::voronoicell_neighbor c;
  	// 	std::vector<std::vector<int> > neighbours(n_positions);
  	// 	std::vector<std::vector<double> > normals(3*n_positions); // inner vector is flattened
  	// 	std::vector<std::vector<double> > face_areas(n_positions);
  	// 	std::vector<double> cell_areas;
  	// 	if( cloop.start() ) do if( con.compute_cell(c,cloop) ){
  	// 		// cloop.pid() gives the index of the cell / particle currently being considered.
  	// 		c.neighbors(  neighbours[cloop.pid()] );
  	// 		c.face_areas( face_areas[cloop.pid()] );
  	// 		c.normals( normals[cloop.pid()] ); // These normals also contain distance, not normalized!
  	// 		cell_areas.push_back( c.surface_area() );
  	// 	} while(cloop.inc());
    //
  	// 	// for (size_t i = 0; i < neighbours[0].size(); i++) {
  	// 	// 	std::cout << neighbours[0][i] << '\n';
  	// 	// }
  	// 	// std::cout << '\n';
    //
    //
  	// 	unsigned int lmax = 6; // TODO replace this by summing over all _requested_ l
  	// 	for(size_t i = 0; i < n_positions; i++){
  	// 		for(size_t l = 0; l < lmax; l++){ q[i][l] = 0.0; } //TEMP
  	// 	}
    //
  	// 	// With neighbours and weights for the MSM, we can compute them
  	// 	for(size_t i = 0; i < n_positions; i++){
  	// 		// Loop over all l that we want to know
  	// 		for(unsigned int l = 0; l < lmax; l++){
  	// 			double factor = 4 * M_PI / (2*l + 1);
  	// 			// Loop over the neighbours / facets
  	// 			for(size_t j = 0; j < neighbours[i].size(); j++){
  	// 				// Calculate the weight factor from the area
  	// 				double area_weight = face_areas[i][j] / cell_areas[i];
  	// 				// Calculate the spherical harmonic
  	// 				for(int m = -((int) l); m <= ((int) l); m++){
  	// 					double dx = positions[i][0] - positions[ neighbours[i][j] ][0];
  	// 					double dy = positions[i][1] - positions[ neighbours[i][j] ][1];
  	// 					// theta (polar) [0,pi], phi (azimuthal) [-pi,pi]
  	// 					double costh = positions[i][2] - positions[ neighbours[i][j] ][2]; // cos(theta) = dz
  	// 					double phi = atan2(dy,dx);
  	// 					std::complex<double> Ylm = spherical_harmonic(l, m, costh, phi); // TODO optimize angles
  	// 					q[i][l] += pow(std::abs(area_weight * Ylm), 2);
  	// 				}
  	// 			}
  	// 			q[i][l] = sqrt(factor * q[i][l]);
  	// 		}
  	// 	}
  	// }
}
