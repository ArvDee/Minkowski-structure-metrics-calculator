#include "process.h"

// Initialize from a unit cell / particle configuration file
void SnapshotProcessor::load_snapshot(const std::string file_name){
	// Open initial config file for reading and feed all data to an array of strings
  std::vector<std::string> text;
  std::string line;
  std::ifstream file(file_name);
  if( file.is_open() ){
		while( getline(file, line) ){
			text.push_back(line);
		}
    file.close();
		if(text.size() < 3){
			std::cout << "Bad initial file, not enough data or wrong format." << '\n';
		}
  } else std::cout << "Unable to open file.\n";

	// Now convert this data to the correct format
	int N = strtol(text[0].c_str(), NULL, 10); // # of particles
  // columns are lattice vecs
	// sscanf(text[1].c_str(),"%f %f %f %f %f %f %f %f %f",&a1[0],&a1[1],&a1[2],&a2[0],&a2[1],&a2[2],&a3[0],&a3[1],&a3[2]); // col-major
	sscanf(text[1].c_str(),"%f %f %f %f %f %f %f %f %f",&a1[0],&a2[0],&a3[0],&a1[1],&a2[1],&a3[1],&a1[2],&a2[2],&a3[2]); // row-major
  positions.reserve(N);
	for(int i = 2; i < N+2; i++){
		// Set the particle's position and orientation from the file
		std::vector<float> pos(3);
		sscanf(text[i].c_str(),"%f %f %f",&pos[0],&pos[1],&pos[2]);
		positions.push_back(pos);
	}
}

// Calculates the bond order parameters ql and wl
void SnapshotProcessor::calculate_order_parameters(size_t max_l){
  size_t N = positions.size();
  // Calculate the bond order parameters
  MSM::MinkowskiStructureCalculator msm;
  msm.load_configuration(positions, a1, a2, a3);
  q = std::vector<std::vector<float>>(max_l+1, std::vector<float>(N));
  w = std::vector<std::vector<float>>(max_l+1, std::vector<float>(N));
  for(size_t l = 0; l <= max_l; l++){
    q[l] = msm.ql_av_all(l);
    w[l] = msm.wl_av_all(l);
  }

  // test
  // float dot = msm.bond_crystallinity(0, 1, 6);
  // std::cout << dot << '\n';
  // std::vector<unsigned int> all_l(max_l-1); // leave out l=0 and l=1
  // for(unsigned int i = 0, l = 2; l < max_l; i++, l++){
  //   all_l[i] = l;
  // }
  // float dot_av = msm.bond_crystallinity_lav(0, 1, all_l);
  // std::cout << dot_av << '\n';
}

// Short test function to check whether a certain file already exists
bool SnapshotProcessor::file_exists(const std::string& name)const{
    std::ifstream f(name.c_str());
    return f.good();
}

// Creates new files to save the q and w bond order parameters to
void SnapshotProcessor::save_qw_files(
  const std::string target_dir,
  const std::string q_file_name = "q.txt",
  const std::string w_file_name = "w.txt"
)const{
  size_t N = positions.size();
  size_t max_l = q.size()-1;
  unsigned int idx = 2;
  std::string q_file_name_ = target_dir + q_file_name;
  std::string w_file_name_ = target_dir + w_file_name;
  // Check if the file name we want already exists
  while( file_exists(q_file_name) || file_exists(w_file_name) ){
    std::string extra_name = "_" + std::to_string(idx);
    q_file_name_ = target_dir + q_file_name + extra_name;
    w_file_name_ = target_dir + w_file_name + extra_name;
    idx++;
  }
  std::ofstream file;
  // Save bond order parameters q of only the particles in the unit cell to a file
  file.open(q_file_name, std::ios::out | std::ios::out);
  for(size_t i = 0; i < N; i++){
    for(size_t l = 0; l <= max_l; l++){file << q[l][i] << " ";}
    file << '\n';
  }
  std::cout << "Saved ql data to " << q_file_name << '\n';
  file.close();
  file.open(w_file_name, std::ios::out | std::ios::out);
  for(size_t i = 0; i < N; i++){
    for(size_t l = 0; l <= max_l; l++){file << w[l][i] << " ";}
    file << '\n';
  }
  std::cout << "Saved wl data to " << w_file_name << '\n';
  file.close();
}
