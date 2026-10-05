#pragma once

// Runs a CaseConfig: conduction first (if configured), then the structural model
// loaded by the resulting temperature, then writes HDF5/XDMF results.

#include <bladefem/case_config.hpp>
#include <bladefem/mesh.hpp>

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace bladefem {

struct CaseReport {
  std::size_t nodes = 0, triangles = 0;

  bool thermal = false;
  double T_min = 0.0, T_max = 0.0;
  Vec2 T_max_at = Vec2::Zero();
  std::vector<std::pair<std::string, double>> boundary_heat;  // [W/m], into the body

  bool structural = false;
  double max_von_mises = 0.0;
  Vec2 max_von_mises_at = Vec2::Zero();
  double min_stress_zz = 0.0, max_stress_zz = 0.0;
  double max_displacement = 0.0;
  Eigen::Vector3d axial_strain = Eigen::Vector3d::Zero();
  double max_reaction = 0.0;  // largest support reaction [N/m]

  std::vector<std::string> files;
};

CaseReport run_case(const CaseConfig& config, const Mesh& mesh);
void print_report(std::ostream& out, const CaseConfig& config, const CaseReport& report);

}  // namespace bladefem
