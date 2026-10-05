#include <bladefem/elastic2d.hpp>
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

Mesh distorted_square(int n, double amount) {
  Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, n, n);
  const double h = 1.0 / n;
  for (std::size_t i = 0; i < m.node_count(); ++i) {
    Vec2& p = m.nodes[i];
    if (p.x() < 1e-12 || p.x() > 1.0 - 1e-12 || p.y() < 1e-12 || p.y() > 1.0 - 1e-12) continue;
    p += amount * h * Vec2(std::sin(12.9898 * static_cast<double>(i) + 78.233), std::cos(39.346 * static_cast<double>(i) + 11.135));
  }
  return m;
}

// Relative L2 error of the element stresses (constant per CST) against an exact field.
template <class Exact>
double stress_l2_relative(const Mesh& m, const ElasticSolution& s, Exact exact) {
  double num = 0.0, den = 0.0;
  for (std::size_t t = 0; t < m.triangle_count(); ++t) {
    const Eigen::Vector3d e = exact(m.centroid(t));
    num += (s.stress[t] - e).squaredNorm() * m.area(t);
    den += e.squaredNorm() * m.area(t);
  }
  return std::sqrt(num / den);
}

}  // namespace

SCENARIO("Patch test: a linear displacement field gives exact uniform stress", "[elastic]") {
  GIVEN("a distorted mesh with u = (1e-3 x + 2e-4 y, -5e-4 x + 3e-4 y) on its boundary") {
    const Mesh m = distorted_square(6, 0.25);
    const auto ux = [](const Vec2& p) { return 1e-3 * p.x() + 2e-4 * p.y(); };
    const auto uy = [](const Vec2& p) { return -5e-4 * p.x() + 3e-4 * p.y(); };

    ElasticProblem prob;
    prob.material = {200e9, 0.3, 0.0};
    for (const char* side : {"bottom", "right", "top", "left"}) {
      prob.displacement.push_back({side, 0, ux});
      prob.displacement.push_back({side, 1, uy});
    }

    WHEN("it is solved in plane stress") {
      const ElasticSolution s = solve_elastic(m, prob);

      THEN("interior displacements and all element stresses are exact") {
        for (std::size_t i = 0; i < m.node_count(); ++i) {
          const Vec2 u = s.node_displacement(static_cast<int>(i));
          CHECK_THAT(u.x(), WithinAbs(ux(m.nodes[i]), 1e-15));
          CHECK_THAT(u.y(), WithinAbs(uy(m.nodes[i]), 1e-15));
        }
        const double E = 200e9, nu = 0.3, c = E / (1 - nu * nu);
        const double ex = 1e-3, ey = 3e-4, gxy = 2e-4 - 5e-4;
        for (const auto& sig : s.stress) {
          CHECK_THAT(sig[0], WithinRel(c * (ex + nu * ey), 1e-9));
          CHECK_THAT(sig[1], WithinRel(c * (ey + nu * ex), 1e-9));
          CHECK_THAT(sig[2], WithinRel(E / (2 * (1 + nu)) * gxy, 1e-9));
        }
      }
    }
  }
}

SCENARIO("A bar in uniaxial tension", "[elastic]") {
  GIVEN("a 2 x 0.5 plate, roller on the left, 100 MPa traction on the right") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 2.0, 0.5, 8, 4);
    ElasticProblem prob;
    prob.material = {70e9, 0.33, 0.0};
    prob.displacement.push_back({"left", 0, [](const Vec2&) { return 0.0; }});
    prob.supports.push_back({Vec2(0.0, 0.0), false, true});
    prob.traction.push_back({"right", [](const Vec2&) { return Vec2(100e6, 0.0); }});

    WHEN("it is solved") {
      const ElasticSolution s = solve_elastic(m, prob);

      THEN("the stress is uniform and the plate stretches and narrows by Hooke's law") {
        for (const auto& sig : s.stress) {
          CHECK_THAT(sig[0], WithinRel(100e6, 1e-9));
          CHECK_THAT(sig[1], WithinAbs(0.0, 1e-3));
        }
        const int corner = m.nearest_node(Vec2(2.0, 0.5));
        CHECK_THAT(s.node_displacement(corner).x(), WithinRel(100e6 / 70e9 * 2.0, 1e-9));
        CHECK_THAT(s.node_displacement(corner).y(), WithinRel(-0.33 * 100e6 / 70e9 * 0.5, 1e-9));
      }
      THEN("the support reaction balances the applied load") {
        CHECK_THAT(s.total_reaction().x(), WithinRel(-100e6 * 0.5, 1e-9));
        CHECK_THAT(s.total_reaction().y(), WithinAbs(0.0, 1e-3));
      }
      THEN("von Mises equals the uniaxial stress") {
        for (double v : s.von_mises) CHECK_THAT(v, WithinRel(100e6, 1e-9));
      }
    }
  }
}

SCENARIO("Timoshenko cantilever: convergence to the elasticity solution", "[elastic][convergence]") {
  GIVEN("L = 48, D = 12, E = 3e7, nu = 0.3, P = 1000 with exact end conditions") {
    const verification::TimoshenkoBeam beam{48.0, 12.0, 3.0e7, 0.3, 1000.0};

    struct Level { double h, disp_l2, tip_error; };
    auto run = [&](int nx) {
      const Mesh m = meshgen::rectangle(0.0, -beam.D / 2, beam.L, beam.D / 2, nx, nx / 4);
      ElasticProblem prob;
      prob.material = {beam.E, beam.nu, 0.0};
      prob.displacement.push_back({"left", 0, [&](const Vec2& p) { return beam.displacement(p).x(); }});
      prob.displacement.push_back({"left", 1, [&](const Vec2& p) { return beam.displacement(p).y(); }});
      prob.traction.push_back({"right", [&](const Vec2& p) {
                                 const auto sig = beam.stress(p);
                                 return Vec2(sig[0], sig[2]);
                               }});
      const ElasticSolution s = solve_elastic(m, prob);

      Eigen::VectorXd ux(m.node_count()), uy(m.node_count());
      for (std::size_t i = 0; i < m.node_count(); ++i) {
        ux[static_cast<Eigen::Index>(i)] = s.displacement[static_cast<Eigen::Index>(2 * i)];
        uy[static_cast<Eigen::Index>(i)] = s.displacement[static_cast<Eigen::Index>(2 * i + 1)];
      }
      const double ex = verification::l2_error(m, ux, [&](const Vec2& p) { return beam.displacement(p).x(); });
      const double ey = verification::l2_error(m, uy, [&](const Vec2& p) { return beam.displacement(p).y(); });
      const double norm = std::hypot(verification::l2_norm(m, [&](const Vec2& p) { return beam.displacement(p).x(); }),
                                     verification::l2_norm(m, [&](const Vec2& p) { return beam.displacement(p).y(); }));
      const double tip = s.node_displacement(m.nearest_node(Vec2(beam.L, 0.0))).y();
      return Level{m.max_edge_length(), std::hypot(ex, ey) / norm, std::abs(tip / beam.tip_deflection() - 1.0)};
    };

    WHEN("the mesh is refined from 16 x 4 to 128 x 32") {
      const Level a = run(16), b = run(32), c = run(64), d = run(128);

      THEN("the displacement error converges at second order") {
        CHECK(verification::observed_rate(c.disp_l2, d.disp_l2, c.h, d.h) == Catch::Approx(2.0).margin(0.1));
      }
      THEN("the tip deflection approaches the exact value from below (CST is stiff)") {
        CHECK(a.tip_error > b.tip_error);
        CHECK(b.tip_error > c.tip_error);
        CHECK(d.tip_error < 0.01);
      }
    }
  }
}

SCENARIO("Kirsch plate: stress concentration at a circular hole", "[elastic][convergence]") {
  GIVEN("a quarter plate a = 1, L = 5 loaded by the exact Kirsch tractions") {
    const verification::KirschPlate kirsch{1.0, 1.0};

    auto run = [&](int n) {
      const Mesh m = meshgen::quarter_plate_with_hole(1.0, 5.0, n, 3 * n / 2, 2.0);
      ElasticProblem prob;
      prob.material = {1000.0, 0.3, 0.0};
      prob.displacement.push_back({"symmetry_x", 0, [](const Vec2&) { return 0.0; }});
      prob.displacement.push_back({"symmetry_y", 1, [](const Vec2&) { return 0.0; }});
      prob.traction.push_back({"right", [&](const Vec2& p) {
                                 const auto sig = kirsch.stress(p);
                                 return Vec2(sig[0], sig[2]);
                               }});
      prob.traction.push_back({"top", [&](const Vec2& p) {
                                 const auto sig = kirsch.stress(p);
                                 return Vec2(sig[2], sig[1]);
                               }});
      const ElasticSolution s = solve_elastic(m, prob);

      std::vector<double> sxx(m.triangle_count());
      for (std::size_t t = 0; t < m.triangle_count(); ++t) sxx[t] = s.stress[t][0];
      const double scf = nodal_average(m, sxx)[static_cast<std::size_t>(m.nearest_node(Vec2(0.0, 1.0)))];
      return std::pair{scf, stress_l2_relative(m, s, [&](const Vec2& p) { return kirsch.stress(p); })};
    };

    WHEN("the mesh is refined") {
      const auto [scf8, e8] = run(8);
      const auto [scf16, e16] = run(16);
      const auto [scf32, e32] = run(32);

      THEN("the stress error halves with each refinement (first order for CST)") {
        CHECK(std::log2(e16 / e32) == Catch::Approx(1.0).margin(0.15));
        CHECK(e8 > e16);
      }
      THEN("the concentration factor at the top of the hole approaches 3") {
        CHECK(std::abs(scf32 - 3.0) < std::abs(scf8 - 3.0));
        CHECK_THAT(scf32, WithinRel(3.0, 0.05));
      }
    }
  }
}

SCENARIO("A structure without enough supports", "[elastic]") {
  GIVEN("a plate held only in x on its left edge") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 3, 3);
    ElasticProblem prob;
    prob.material = {1e9, 0.3, 0.0};
    prob.displacement.push_back({"left", 0, [](const Vec2&) { return 0.0; }});
    THEN("the remaining rigid-body translation is reported instead of a bogus answer") {
      CHECK_THROWS_AS(solve_elastic(m, prob), std::runtime_error);
    }
  }
  GIVEN("an invalid Poisson's ratio") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 1, 1);
    ElasticProblem prob;
    prob.material = {1e9, 0.5, 0.0};
    THEN("solving is refused") { CHECK_THROWS_AS(solve_elastic(m, prob), std::invalid_argument); }
  }
}
