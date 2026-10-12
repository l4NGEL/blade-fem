#pragma once

// Steady 2D heat conduction on a linear-triangle mesh:
//   -div( k grad T ) = Q            in the domain
//    T = T_D                        on Dirichlet boundaries
//   -k dT/dn = h (T - T_inf)        on convection boundaries
//   -k dT/dn = -q_in                on flux boundaries (q_in > 0 heats the body)

#include <bladefem/mesh.hpp>

#include <Eigen/Core>

#include <functional>
#include <string>
#include <vector>

namespace bladefem {

using ScalarField = std::function<double(const Vec2&)>;

struct ThermalProblem {
  struct Dirichlet {
    std::string boundary;
    ScalarField temperature;
  };
  struct Convection {
    std::string boundary;
    double h;      // [W/(m^2 K)]
    double T_inf;  // [K]
  };
  struct HeatFlux {
    std::string boundary;
    double q_in;  // [W/m^2], into the body
  };

  double conductivity = 1.0;  // k [W/(m K)]
  ScalarField source;         // Q [W/m^3]; empty = none
  std::vector<Dirichlet> dirichlet;
  std::vector<Convection> convection;
  std::vector<HeatFlux> flux;
};

struct ThermalSolution {
  Eigen::VectorXd temperature;  // nodal values
  // Net heat entering the body through each named boundary [W per unit depth].
  // For a converged steady solution these sum to minus the internal generation.
  std::vector<std::pair<std::string, double>> boundary_heat;
};

ThermalSolution solve_thermal(const Mesh& mesh, const ThermalProblem& problem);

}  // namespace bladefem
