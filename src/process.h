#ifndef BOOP_SNAPSHOT_H
#define BOOP_SNAPSHOT_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include "MinkowskiStructureMetrics.h"

class SnapshotProcessor{
private:
public:
  std::vector<float> box = std::vector<float>(9); // Simulation volume
  std::vector<std::vector<float>> positions;       // Positions of particles in box
  std::vector<std::vector<float>> q;
  std::vector<std::vector<float>> w;
	SnapshotProcessor(void);
	~SnapshotProcessor(void){};

	void load_snapshot(const std::string init_file_name);
  void calculate_order_parameters(size_t max_l);
  void save_qw_files(const std::string input_file_name, const std::string optional_file_string);
  bool file_exists(const std::string& name);
};

#endif // end header guard BOOP_SNAPSHOT_H
