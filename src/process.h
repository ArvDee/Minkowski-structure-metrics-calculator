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
  std::vector<double> a1 = std::vector<double>(3); // Lattice vectors
  std::vector<double> a2 = std::vector<double>(3);
  std::vector<double> a3 = std::vector<double>(3);
  std::vector<std::vector<double>> positions; // Positions of particles in box
  std::vector<std::vector<double>> q;
  std::vector<std::vector<double>> w;
  std::vector<std::vector<double>> q_av;
  std::vector<std::vector<double>> w_av;
  std::vector<std::vector<double>> q_dot;
  std::vector<std::vector<double>> q_dot_av;
  SnapshotProcessor(void){};
  virtual ~SnapshotProcessor(void){};

  void load_snapshot(const std::string file_name);

  void calculate_order_parameters(size_t max_l);

  bool file_exists(const std::string& name)const;

  void save_boops(
    const std::string target_dir,
    const std::string file_name,
    std::vector<std::vector<double>>& boop_vector
  )const;
};

#endif // end header guard BOOP_SNAPSHOT_H
