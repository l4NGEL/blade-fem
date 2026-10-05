#include <bladefem/thermal2d.hpp>
#include <bladefem/verification.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

using namespace bladefem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Moves interior nodes by a deterministic pseudo-random offset of up to `amount`
// cell sizes, so that patch tests do not rely on a regular grid.
Mesh distorted_square(int n, double amount) {
  Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, n, n);
  const double h = 1.0 / n;
  for (std::size_t i = 0; i < m.node_count(); ++i) {
    Vec2& p = m.nodes[i];
    if (p.x() < 1e-12 || p.x() > 1.0 - 1e-12 || p.y() < 1e-12 || p.y() > 1.0 - 1e-12) continue;
    const double s = std::sin(12.9898 * static_cast<double>(i) + 78.233);
    const double c = std::cos(39.346 * static_cast<double>(i) + 11.135);
    p += amount * h * Vec2(s, c);
  }
  return m;
}

ThermalProblem all_dirichlet(ScalarField T) {
  ThermalProblem p;
  for (const char* side : {"bottom", "right", "top", "left"}) p.dirichlet.push_back({side, T});
  return p;
}

}  // namespace

SCENARIO("Patch test: a linear temperature field is reproduced exactly", "[thermal]") {
  GIVEN("a distorted mesh with a linear field prescribed on its boundary") {
    const Mesh m = distorted_square(6, 0.25);
    const auto exact = [](const Vec2& p) { return 300.0 + 40.0 * p.x() - 25.0 * p.y(); };
    ThermalProblem prob = all_dirichlet(exact);
    prob.conductivity = 7.5;

    WHEN("it is solved") {
      const ThermalSolution s = solve_thermal(m, prob);
      THEN("every interior node matches the linear field to round-off") {
        for (std::size_t i = 0; i < m.node_count(); ++i) {
          CHECK_THAT(s.temperature[static_cast<Eigen::Index>(i)], WithinAbs(exact(m.nodes[i]), 1e-9));
        }
      }
    }
  }
}

SCENARIO("Manufactured solution: convergence rates of linear triangles", "[thermal][convergence]") {
  GIVEN("T = sin(pi x) sin(pi y) on the unit square with the matching source") {
    const double k = 4.0;
    const auto exact = [](const Vec2& p) { return std::sin(kPi * p.x()) * std::sin(kPi * p.y()); };
    const auto grad = [](const Vec2& p) {
      return Vec2(kPi * std::cos(kPi * p.x()) * std::sin(kPi * p.y()),
                  kPi * std::sin(kPi * p.x()) * std::cos(kPi * p.y()));
    };

    struct Level { double h, l2, h1; };
    auto run = [&](int n) {
      const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, n, n);
      ThermalProblem prob = all_dirichlet([](const Vec2&) { return 0.0; });
      prob.conductivity = k;
      prob.source = [&](const Vec2& p) { return 2.0 * kPi * kPi * k * exact(p); };
      const ThermalSolution s = solve_thermal(m, prob);
      return Level{m.max_edge_length(), verification::l2_error(m, s.temperature, exact),
                   verification::h1_seminorm_error(m, s.temperature, grad)};
    };

    WHEN("the mesh is refined from 8 to 64 cells per side") {
      const Level a = run(8), b = run(16), c = run(32), d = run(64);

      THEN("the L2 error converges at order 2 and the H1 error at order 1") {
        CHECK(verification::observed_rate(c.l2, d.l2, c.h, d.h) == Catch::Approx(2.0).margin(0.05));
        CHECK(verification::observed_rate(c.h1, d.h1, c.h, d.h) == Catch::Approx(1.0).margin(0.05));
        CHECK(a.l2 > b.l2);
        CHECK(b.l2 > c.l2);
      }
    }
  }
}

SCENARIO("Radial conduction through a ring cooled on the outside", "[thermal]") {
  GIVEN("T = 900 K on r = 0.02 m and convection h = 400, T_inf = 300 K on r = 0.05 m") {
    const double ri = 0.02, ro = 0.05, k = 15.0, h = 400.0, Tin = 900.0, Tinf = 300.0;
    // T(r) = T_i + C ln(r / r_i),  -k T'(r_o) = h (T(r_o) - T_inf)
    const double C = h * (Tinf - Tin) / (k / ro + h * std::log(ro / ri));
    const auto exact = [&](const Vec2& p) { return Tin + C * std::log(p.norm() / ri); };

    ThermalProblem prob;
    prob.conductivity = k;
    prob.dirichlet.push_back({"inner", [&](const Vec2&) { return Tin; }});
    prob.convection.push_back({"outer", h, Tinf});

    WHEN("it is solved on a 128 x 16 ring mesh") {
      const Mesh m = meshgen::annulus(ri, ro, 128, 16);
      const ThermalSolution s = solve_thermal(m, prob);

      THEN("the outer surface temperature matches the analytical value") {
        for (int n : m.boundary_nodes("outer")) {
          CHECK_THAT(s.temperature[n], WithinRel(exact(m.nodes[static_cast<std::size_t>(n)]), 2e-3));
        }
      }
      THEN("the heat entering at the inner wall leaves through the outer one") {
        double in = 0.0, out = 0.0;
        for (const auto& [name, q] : s.boundary_heat) {
          if (name == "inner") in = q;
          if (name == "outer") out = q;
        }
        CHECK(in > 0.0);
        CHECK_THAT(in + out, WithinAbs(0.0, 1e-8 * in));
        const double analytical = -k * C / ri * 2.0 * kPi * ri;  // per unit depth
        CHECK_THAT(in, WithinRel(analytical, 1e-2));
      }
    }
  }
}

SCENARIO("A slab heated by a prescribed flux", "[thermal]") {
  GIVEN("q = 2 kW/m^2 entering on the left, T = 300 K on the right, insulated top and bottom") {
    const Mesh m = distorted_square(5, 0.2);
    ThermalProblem prob;
    prob.conductivity = 10.0;
    prob.flux.push_back({"left", 2000.0});
    prob.dirichlet.push_back({"right", [](const Vec2&) { return 300.0; }});

    WHEN("it is solved") {
      const ThermalSolution s = solve_thermal(m, prob);
      THEN("the temperature falls linearly by q L / k") {
        for (std::size_t i = 0; i < m.node_count(); ++i) {
          CHECK_THAT(s.temperature[static_cast<Eigen::Index>(i)],
                     WithinAbs(300.0 + 200.0 * (1.0 - m.nodes[i].x()), 1e-9));
        }
      }
    }
  }
}

SCENARIO("Ill-posed thermal problems are rejected", "[thermal]") {
  const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 2, 2);
  GIVEN("a boundary condition on a boundary that does not exist") {
    ThermalProblem prob;
    prob.convection.push_back({"inlet", 10.0, 300.0});
    THEN("solving throws") { CHECK_THROWS_AS(solve_thermal(m, prob), std::out_of_range); }
  }
  GIVEN("only flux conditions, so the temperature level is free") {
    ThermalProblem prob;
    prob.flux.push_back({"left", 1.0});
    THEN("solving throws") { CHECK_THROWS_AS(solve_thermal(m, prob), std::runtime_error); }
  }
}
