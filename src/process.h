#ifndef BOOP_SNAPSHOT_H
#define BOOP_SNAPSHOT_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <eigen3/Eigen/Geometry>
#include "MinkowskiStructureMetrics.h"

class SnapshotProcessor{
private:
public:
  std::vector<float> a1 = std::vector<float>(3); // Lattice vectors
  std::vector<float> a2 = std::vector<float>(3);
  std::vector<float> a3 = std::vector<float>(3);
  std::vector<std::vector<float>> positions; // Positions of particles in box
  std::vector<std::vector<float>> q;
  std::vector<std::vector<float>> w;
	SnapshotProcessor(void){};
	virtual ~SnapshotProcessor(void){};

	void load_snapshot(const std::string file_name);
  void calculate_order_parameters(size_t max_l);
  bool file_exists(const std::string& name)const;
  void save_qw_files(const std::string target_dir, const std::string q_file_name, const std::string w_file_name)const;
};

#endif // end header guard BOOP_SNAPSHOT_H
