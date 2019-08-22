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
  q_av = std::vector<std::vector<float>>(max_l+1, std::vector<float>(N));
  w_av = std::vector<std::vector<float>>(max_l+1, std::vector<float>(N));
  for(size_t l = 0; l <= max_l; l++){
    q[l] = msm.ql_all(l);
    w[l] = msm.wl_all(l);
    q_av[l] = msm.ql_av_all(l);
    w_av[l] = msm.wl_av_all(l);
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
    return f.good(); // f is closed automatically on function return
}

// Creates new files to save the bond order parameters to
void SnapshotProcessor::save_boops(
  const std::string target_dir,
  const std::string file_name,
  std::vector<std::vector<float>>& boop_vector
)const{
  size_t N = positions.size();
  size_t max_l = boop_vector.size()-1;
  std::string extension = ".txt";
  unsigned int idx = 0;
  // Create the full path
  std::string full_path = target_dir + file_name + '_' + std::to_string(idx) + extension;
  // Check if the file name we want already exists
  while( file_exists(full_path) ){
    idx++;
    full_path = target_dir + file_name + '_' + std::to_string(idx) + extension;
  }
  std::ofstream file;
  // Write data
  file.open(full_path, std::ios::out | std::ios::out);
  for(size_t i = 0; i < N; i++){
    for(size_t l = 0; l <= max_l; l++){file << boop_vector[l][i] << " ";}
    file << '\n';
  }
  std::cout << "Saved boop data to " << full_path << '\n';
  file.close();
}
