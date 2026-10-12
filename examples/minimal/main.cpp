// Smallest useful program against the library: steady conduction through a
// ring heated on the inside and cooled by convection on the outside, compared
// with the closed-form solution.

#include <bladefem/mesh.hpp>
#include <bladefem/thermal2d.hpp>
#include <bladefem/version.hpp>

#include <cmath>
#include <cstdio>

int main() {
  using namespace bladefem;

  const double r_in = 0.02, r_out = 0.05, k = 15.0, h = 400.0, T_in = 900.0, T_inf = 300.0;
  const Mesh mesh = meshgen::annulus(r_in, r_out, 128, 16);

  ThermalProblem problem;
  problem.conductivity = k;
  problem.dirichlet.push_back({"inner", [&](const Vec2&) { return T_in; }});
  problem.convection.push_back({"outer", h, T_inf});
  const ThermalSolution solution = solve_thermal(mesh, problem);

  const double C = h * (T_inf - T_in) / (k / r_out + h * std::log(r_out / r_in));
  const double exact = T_in + C * std::log(r_out / r_in);
  const int probe = mesh.nearest_node(Vec2(r_out, 0.0));

  std::printf("blade-fem %s\n", kVersion);
  std::printf("outer wall: %.2f K (FE) vs %.2f K (exact)\n", solution.temperature[probe], exact);
  return std::abs(solution.temperature[probe] - exact) < 1.0 ? 0 : 1;
}
