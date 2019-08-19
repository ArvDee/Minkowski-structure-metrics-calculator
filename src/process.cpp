#include "process.h"


/* ------------- Constructor -------------
 *
 */
SnapshotProcessor::SnapshotProcessor(void){
}

/* --------- Initialize from a unit cell / particle configuration file ---------
 *
 */
void SnapshotProcessor::load_snapshot(const std::string init_file_name){
	// Open initial config file for reading and feed all data to an array of strings
	std::vector<std::string> text;
	std::string line;
  std::ifstream init_file(init_file_name);
  if( init_file.is_open() ){
		while( getline(init_file, line) ){
			text.push_back(line);
		}
    init_file.close();
		if(text.size() < 3){
			std::cout << "Bad initial file, not enough data or wrong format." << '\n';
		}
  } else std::cout << "Unable to open file.\n";

	// Now convert this data to the correct format
	int N = strtol(text[0].c_str(), NULL, 10); // # of particles
	sscanf(text[1].c_str(),"%f %f %f %f %f %f %f %f %f",&box[0],&box[1],&box[2],&box[3],&box[4],&box[5],&box[6],&box[7],&box[8]); // col-major
	// sscanf(text[1].c_str(),"%f %f %f %f %f %f %f %f %f",&box[0],&box[3],&box[6], &box[1],&box[4],&box[7], &box[2],&box[5],&box[8]); // row-major
  // for(auto el : box){
  //   std::cout << el << ' ';
  // }
  // std::cout << '\n';
  positions.reserve(N);
	for(int i = 0; i < N; i++){
		// Set the particle's position and orientation from the file
		std::vector<float> pos(3);
    float angle, axis[3];
		sscanf(text[i+2].c_str(),"%f %f %f %f %f %f %f",&pos[0],&pos[1],&pos[2],&angle,&axis[1],&axis[2],&axis[3]);
		positions.push_back(pos);
	}

	// for(int i = 0; i < N; i++){
  //   std::cout << positions[i][0] <<' '<< positions[i][1] <<' '<< positions[i][2] << '\n';
  // }
	// for(int i = 0; i < 9; i++){
  //   std::cout << box[i] << ' ';
  // }
  // std::cout << '\n';
}

/* --------- Calculates the bond order parameters ql and wl ---------
 * Because Voro++ only works with upper triangular box matrices, there are
 * some extra steps.
 */
void SnapshotProcessor::calculate_order_parameters_generalBox(const size_t max_l){
  // For all particles, generate copies inside a cube
  // size_t copies = 5 - std::min( floor(0.5*pow(positions.size(),1.0/3.0)), 4.0);
  size_t copies = 1;
  std::vector<std::vector<float>> positions_with_copies = positions;
  std::vector<float> pos(3), rel_pos(3);
  for(size_t idx = 0; idx < positions.size(); idx++){
    // Then make new positions by copying the unit cell
    for(int i = -((int) copies); i <= ((int) copies); i++){
      for(int j = -((int) copies); j <= ((int) copies); j++){
        for(int k = -((int) copies); k <= ((int) copies); k++){
          if(i==0 && j==0 && k==0){ continue; } // Skip position of particle itself
          rel_pos[0] = float(i);
          rel_pos[1] = float(j);
          rel_pos[2] = float(k);
          // Create new positions at integer multiples of the box vectors (columns)
          pos[0] = positions[idx][0] + box[0] * rel_pos[0] + box[1] * rel_pos[1] + box[2] * rel_pos[2];
          pos[1] = positions[idx][1] + box[3] * rel_pos[0] + box[4] * rel_pos[1] + box[5] * rel_pos[2];
          pos[2] = positions[idx][2] + box[6] * rel_pos[0] + box[7] * rel_pos[1] + box[8] * rel_pos[2];
          positions_with_copies.push_back(pos);
        }
      }
    }
  }
  std::vector<float> minpos = positions[0];
  std::vector<float> maxpos = positions[0];
  // Find the position that will mark the corner of the box with copies
  for(size_t i = 0; i < positions_with_copies.size(); i++){
    for(size_t d = 0; d < 3; d++){
      if(positions_with_copies[i][d] < minpos[d]){ minpos[d] = positions_with_copies[i][d]; }
      if(positions_with_copies[i][d] > maxpos[d]){ maxpos[d] = positions_with_copies[i][d]; }
    }
  }
  // Offset to put all positions in positive octant
  for(size_t i = 0; i < positions_with_copies.size(); i++){
     for(size_t d = 0; d < 3; d++){
       positions_with_copies[i][d] -= minpos[d];
     }
  }
  // Create a cubic box that can contain all these neighbours
  std::vector<float> boxm(9, 0.0);
  for(size_t i = 0; i < 3; i++){
    boxm[4*i] = maxpos[i]-minpos[i]+1E-10; // offset to prevent edge cases
  }
  // Calculate the bond order parameters for only the original particles
  q = std::vector<std::vector<float>>(positions.size(), std::vector<float>(max_l+1));
  w = std::vector<std::vector<float>>(positions.size(), std::vector<float>(max_l+1));
  // MSM::msm(positions_with_copies, boxm, q, w);

  MSM::MinkowskiStructureCalculator msm;
  msm.compute(positions, box, q, w);
}

/* --------- Calculates the bond order parameters ql and wl ---------
 *
 */
void SnapshotProcessor::calculate_order_parameters_upperTriangularBox(const size_t max_l){
  size_t N = positions.size();
  // Calculate the bond order parameters
  q = std::vector<std::vector<float>>(N, std::vector<float>(max_l+1));
  w = std::vector<std::vector<float>>(N, std::vector<float>(max_l+1));
  MSM::MinkowskiStructureCalculator msm;
  msm.compute(positions, box, q, w);
}

// Short test function to check whether a certain file already exists
bool SnapshotProcessor::file_exists(const std::string& name){
    std::ifstream f(name.c_str());
    return f.good();
}

/* --------- Creates new files to save the q and w bond order parameters to ---------
 *
 */
void SnapshotProcessor::save_qw_files(const std::string target_dir, const std::string optional_file_string){
  size_t N = positions.size();
  size_t max_l = q[0].size()-1;
  unsigned int idx = 2;
  std::string q_file_name = target_dir + "q" + optional_file_string + ".txt";
  std::string w_file_name = target_dir + "w" + optional_file_string + ".txt";
  // Check if the file name we want already exists
  while( file_exists(q_file_name) || file_exists(w_file_name) ){
    std::string extra_name = "_" + std::to_string(idx);
    q_file_name = target_dir + "q" + optional_file_string + extra_name + ".txt";
    w_file_name = target_dir + "w" + optional_file_string + extra_name + ".txt";
    idx++;
  }
  std::ofstream file;
  // Save bond order parameters q of only the particles in the unit cell to a file
  file.open(q_file_name, std::ios::out | std::ios::out);
  for(size_t i = 0; i < N; i++){
    for(size_t j = 0; j <= max_l; j++){file << q[i][j] << " ";}
    file << '\n';
  }
  std::cout << "Saved ql data to " << q_file_name << '\n';
  file.close();
  file.open(w_file_name, std::ios::out | std::ios::out);
  for(size_t i = 0; i < N; i++){
    for(size_t j = 0; j <= max_l; j++){file << w[i][j] << " ";}
    file << '\n';
  }
  std::cout << "Saved wl data to " << w_file_name << '\n';
  file.close();
}
