#include <bladefem/fem1d.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

using namespace bladefem::fem1d;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double kPi = 3.14159265358979323846;

// L2 norm of (T_h - T_exact) over [0, L], sampled densely enough that the
// sampling error is far below the discretisation error being measured.
template <class Exact>
double l2_error(const BarSolution& s, double L, Exact exact) {
  const int samples = 20000;
  const double dx = L / samples;
  double sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double x = (i + 0.5) * dx;
    const double e = s.at(x) - exact(x);
    sum += e * e * dx;
  }
  return std::sqrt(sum);
}

}  // namespace

SCENARIO("A bar with a uniform heat source and fixed end temperatures", "[fem1d]") {
  GIVEN("k = 15, Q = 2e4, T(0) = 300 K and T(L) = 400 K on a 0.2 m bar") {
    BarProblem p;
    p.length = 0.2;
    p.elements = 8;
    p.conductivity = 15.0;
    p.source = [](double) { return 2.0e4; };
    p.left = EndCondition::temperature(300.0);
    p.right = EndCondition::temperature(400.0);

    WHEN("it is solved with linear elements") {
      const BarSolution s = solve(p);

      THEN("the nodal temperatures match the parabolic exact solution") {
        const auto exact = [&](double x) {
          return 300.0 + 100.0 * x / p.length + 2.0e4 / (2.0 * p.conductivity) * x * (p.length - x);
        };
        REQUIRE(s.x.size() == 9);
        for (std::size_t i = 0; i < s.x.size(); ++i) {
          CHECK_THAT(s.temperature[i], WithinAbs(exact(s.x[i]), 1e-9));
        }
      }
    }
  }
}

SCENARIO("A bar cooled by convection at one end", "[fem1d]") {
  GIVEN("T(0) = 500 K, no source, and h = 80, T_inf = 300 K at x = L") {
    BarProblem p;
    p.length = 0.5;
    p.elements = 5;
    p.conductivity = 20.0;
    p.left = EndCondition::temperature(500.0);
    p.right = EndCondition::convection(80.0, 300.0);

    WHEN("it is solved") {
      const BarSolution s = solve(p);

      THEN("the linear profile balances conduction against convection") {
        // T = T0 + C x with -k C = h (T(L) - T_inf)
        const double C = 80.0 * (300.0 - 500.0) / (20.0 + 80.0 * 0.5);
        for (std::size_t i = 0; i < s.x.size(); ++i) {
          CHECK_THAT(s.temperature[i], WithinAbs(500.0 + C * s.x[i], 1e-9));
        }
        const double conducted = -20.0 * C;
        const double convected = 80.0 * (s.temperature.back() - 300.0);
        CHECK_THAT(conducted, WithinRel(convected, 1e-12));
      }
    }
  }
}

SCENARIO("A bar heated by a prescribed flux", "[fem1d]") {
  GIVEN("5 kW/m^2 entering at x = 0 and T(L) = 350 K") {
    BarProblem p;
    p.length = 0.1;
    p.elements = 4;
    p.conductivity = 25.0;
    p.left = EndCondition::heat_flux(5000.0);
    p.right = EndCondition::temperature(350.0);

    WHEN("it is solved") {
      const BarSolution s = solve(p);

      THEN("the gradient carries exactly that flux towards the cold end") {
        for (std::size_t i = 0; i < s.x.size(); ++i) {
          CHECK_THAT(s.temperature[i], WithinAbs(350.0 + 5000.0 / 25.0 * (0.1 - s.x[i]), 1e-9));
        }
      }
    }
  }
}

SCENARIO("Mesh refinement on a smooth solution", "[fem1d][convergence]") {
  GIVEN("Q = k (pi/L)^2 sin(pi x / L) with homogeneous end temperatures") {
    const double L = 1.0;
    const auto exact = [&](double x) { return std::sin(kPi * x / L); };

    auto run = [&](int n) {
      BarProblem p;
      p.length = L;
      p.elements = n;
      p.conductivity = 3.0;
      p.source = [&](double x) { return 3.0 * (kPi / L) * (kPi / L) * std::sin(kPi * x / L); };
      p.left = EndCondition::temperature(0.0);
      p.right = EndCondition::temperature(0.0);
      return l2_error(solve(p), L, exact);
    };

    WHEN("the element count is doubled twice") {
      const double e1 = run(10), e2 = run(20), e3 = run(40);

      THEN("the L2 error falls at the theoretical second-order rate") {
        const double r1 = std::log2(e1 / e2);
        const double r2 = std::log2(e2 / e3);
        CHECK(r1 == Catch::Approx(2.0).margin(0.05));
        CHECK(r2 == Catch::Approx(2.0).margin(0.05));
      }
    }
  }
}

SCENARIO("Invalid bar definitions are rejected", "[fem1d]") {
  GIVEN("a bar with no fixed temperature and no convection") {
    BarProblem p;
    p.left = EndCondition::heat_flux(1.0);
    p.right = EndCondition::heat_flux(0.0);
    THEN("solving it throws, because the temperature level is undetermined") {
      CHECK_THROWS_AS(solve(p), std::invalid_argument);
    }
  }
  GIVEN("a bar with zero elements") {
    BarProblem p;
    p.elements = 0;
    p.left = EndCondition::temperature(0.0);
    THEN("solving it throws") { CHECK_THROWS_AS(solve(p), std::invalid_argument); }
  }
}
