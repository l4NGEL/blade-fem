#pragma once

// Linear elasticity with constant-strain triangles (CST), plane stress or plane
// strain, optionally loaded by a nodal temperature field (thermal strain).
// Results are per unit thickness.

#include <bladefem/mesh.hpp>

#include <Eigen/Core>

#include <functional>
#include <string>
#include <vector>

namespace bladefem {

struct ElasticMaterial {
  double E = 0.0;      // Young's modulus [Pa]
  double nu = 0.0;     // Poisson's ratio
  double alpha = 0.0;  // thermal expansion coefficient [1/K]
};

enum class PlaneModel { Stress, Strain };

struct ElasticProblem {
  struct Displacement {
    std::string boundary;
    int component;                              // 0 = x, 1 = y
    std::function<double(const Vec2&)> value;   // [m]
  };
  // Fixes the node nearest to `at`; used to remove rigid-body motion of free bodies.
  struct PointSupport {
    Vec2 at;
    bool fix_x = true;
    bool fix_y = true;
  };
  struct Traction {
    std::string boundary;
    std::function<Vec2(const Vec2&)> traction;  // surface force [Pa]
  };

  ElasticMaterial material;
  PlaneModel model = PlaneModel::Stress;
  std::vector<Displacement> displacement;
  std::vector<PointSupport> supports;
  std::vector<Traction> traction;

  // Nodal temperature [K]; empty = no thermal load. Strain-free at T_ref.
  Eigen::VectorXd temperature;
  double T_ref = 0.0;
};

struct ElasticSolution {
  Eigen::VectorXd displacement;         // 2N, interleaved (ux0, uy0, ux1, ...)
  std::vector<Eigen::Vector3d> stress;  // per triangle: sxx, syy, sxy
  std::vector<double> stress_zz;        // per triangle (non-zero in plane strain only)
  std::vector<double> von_mises;        // per triangle
  Eigen::VectorXd reactions;            // 2N, non-zero at constrained dofs

  Vec2 node_displacement(int n) const { return {displacement[2 * n], displacement[2 * n + 1]}; }
  Vec2 total_reaction() const;
};

ElasticSolution solve_elastic(const Mesh& mesh, const ElasticProblem& problem);

// von Mises stress from in-plane components and the out-of-plane normal stress.
double von_mises(const Eigen::Vector3d& s, double szz = 0.0);

// Area-weighted average of per-triangle values onto the nodes (for plotting and
// for reading boundary values such as stress concentrations).
std::vector<double> nodal_average(const Mesh& mesh, const std::vector<double>& cell_values);

}  // namespace bladefem
