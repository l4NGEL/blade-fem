#pragma once

// Steady 1D heat conduction in a bar of unit cross-section,
//   -d/dx( k dT/dx ) = Q(x),   0 < x < L,
// discretised with linear (2-node) elements.

#include <functional>
#include <vector>

namespace bladefem::fem1d {

// Boundary condition applied at one end of the bar.
struct EndCondition {
  enum class Kind { Temperature, HeatFlux, Convection };

  Kind kind = Kind::HeatFlux;
  double value = 0.0;  // T [K] | heat flux into the bar [W/m^2] | ambient T_inf [K]
  double h = 0.0;      // film coefficient [W/(m^2 K)], Convection only

  static EndCondition temperature(double T);
  static EndCondition heat_flux(double q_in);
  static EndCondition convection(double h, double T_inf);
};

struct BarProblem {
  double length = 1.0;
  int elements = 10;
  double conductivity = 1.0;              // k [W/(m K)]
  std::function<double(double)> source;   // Q(x) [W/m^3]; empty means no source
  EndCondition left = EndCondition::heat_flux(0.0);
  EndCondition right = EndCondition::heat_flux(0.0);
};

struct BarSolution {
  std::vector<double> x;
  std::vector<double> temperature;

  // Piecewise-linear interpolation of the finite element solution.
  double at(double xq) const;
};

BarSolution solve(const BarProblem& problem);

}  // namespace bladefem::fem1d
