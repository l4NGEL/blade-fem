#include <bladefem/elastic2d.hpp>
#include <bladefem/thermal2d.hpp>
#include <bladefem/verification.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace bladefem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double E = 200e9, nu = 0.3, alpha = 13e-6;

Eigen::VectorXd nodal(const Mesh& m, const std::function<double(const Vec2&)>& f) {
  Eigen::VectorXd v(m.node_count());
  for (std::size_t i = 0; i < m.node_count(); ++i) v[static_cast<Eigen::Index>(i)] = f(m.nodes[i]);
  return v;
}

// Statically determinate supports for a body that should expand freely:
// pin at (x0, y0), roller in y at (x1, y0).
void hold_free(ElasticProblem& p, Vec2 pin, Vec2 roller) {
  p.supports.push_back({pin, true, true});
  p.supports.push_back({roller, false, true});
}

double max_abs_stress(const ElasticSolution& s) {
  double m = 0.0;
  for (const auto& sig : s.stress) m = std::max(m, sig.cwiseAbs().maxCoeff());
  return m;
}

}  // namespace

SCENARIO("Uniform heating of a free plate", "[thermal-stress]") {
  GIVEN("a 0.04 x 0.02 m plate heated uniformly by 500 K, held statically determinate") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 0.04, 0.02, 8, 4);
    ElasticProblem prob;
    prob.material = {E, nu, alpha};
    prob.temperature = nodal(m, [](const Vec2&) { return 800.0; });
    prob.T_ref = 300.0;
    hold_free(prob, Vec2(0, 0), Vec2(0.04, 0));

    WHEN("it is solved") {
      const ElasticSolution s = solve_elastic(m, prob);

      THEN("it expands by alpha * dT without any stress") {
        CHECK(max_abs_stress(s) < 1e-6 * E * alpha * 500.0);
        for (std::size_t i = 0; i < m.node_count(); ++i) {
          const Vec2 u = s.node_displacement(static_cast<int>(i));
          CHECK_THAT(u.x(), WithinAbs(alpha * 500.0 * m.nodes[i].x(), 1e-15));
          CHECK_THAT(u.y(), WithinAbs(alpha * 500.0 * m.nodes[i].y(), 1e-15));
        }
      }
      THEN("the supports carry no load") {
        CHECK(s.total_reaction().norm() < 1e-6);
        CHECK(s.reactions.cwiseAbs().maxCoeff() < 1e-6);
      }
    }
  }
}

SCENARIO("Uniform heating of a fully clamped plate", "[thermal-stress]") {
  GIVEN("all edges fixed and dT = 100 K") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 4, 4);
    ElasticProblem prob;
    prob.material = {E, nu, alpha};
    prob.temperature = nodal(m, [](const Vec2&) { return 400.0; });
    prob.T_ref = 300.0;
    for (const char* side : {"bottom", "right", "top", "left"}) {
      prob.displacement.push_back({side, 0, [](const Vec2&) { return 0.0; }});
      prob.displacement.push_back({side, 1, [](const Vec2&) { return 0.0; }});
    }

    WHEN("modelled as plane stress") {
      const ElasticSolution s = solve_elastic(m, prob);
      THEN("the in-plane stress is -E alpha dT / (1 - nu) everywhere") {
        for (const auto& sig : s.stress) {
          CHECK_THAT(sig[0], WithinRel(-E * alpha * 100.0 / (1 - nu), 1e-9));
          CHECK_THAT(sig[1], WithinRel(-E * alpha * 100.0 / (1 - nu), 1e-9));
          CHECK_THAT(sig[2], WithinAbs(0.0, 1e-3));
        }
      }
    }
    WHEN("modelled as plane strain") {
      prob.model = PlaneModel::Strain;
      const ElasticSolution s = solve_elastic(m, prob);
      THEN("all three normal stresses equal -E alpha dT / (1 - 2 nu)") {
        const double expected = -E * alpha * 100.0 / (1 - 2 * nu);
        for (std::size_t t = 0; t < m.triangle_count(); ++t) {
          CHECK_THAT(s.stress[t][0], WithinRel(expected, 1e-9));
          CHECK_THAT(s.stress_zz[t], WithinRel(expected, 1e-9));
          CHECK_THAT(s.von_mises[t], WithinAbs(0.0, 1e-6 * std::abs(expected)));
        }
      }
    }
  }
}

SCENARIO("A linear temperature gradient in a free plate", "[thermal-stress][convergence]") {
  GIVEN("T rising by 400 K across a free square") {
    auto run = [](int n) {
      const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, n, n);
      ElasticProblem prob;
      prob.material = {E, nu, alpha};
      prob.temperature = nodal(m, [](const Vec2& p) { return 300.0 + 400.0 * p.x(); });
      prob.T_ref = 300.0;
      hold_free(prob, Vec2(0, 0), Vec2(1, 0));
      return solve_elastic(m, prob);
    };

    WHEN("solved on successively finer meshes") {
      const ElasticSolution coarse = run(8), fine = run(16), finer = run(32);

      THEN("the self-equilibrated thermal load leaves the supports unloaded") {
        // Nodal thermal forces are of order E alpha dT L ~ 1e9 N; reactions are round-off.
        CHECK(fine.reactions.cwiseAbs().maxCoeff() < 1e-12 * E * alpha * 400.0);
      }
      THEN("the (exactly zero) stress is approached at first order") {
        // A linear T is compatible in plane stress, but the exact displacement is
        // quadratic, so CST shows a small discretisation stress that halves with h.
        const double scale = E * alpha * 400.0;
        CHECK(max_abs_stress(coarse) / scale < 0.1);
        CHECK(std::log2(max_abs_stress(fine) / max_abs_stress(finer)) == Catch::Approx(1.0).margin(0.2));
      }
    }
  }
}

SCENARIO("Coupled analysis: a ring heated inside and cooled outside", "[thermal-stress][coupled][convergence]") {
  GIVEN("T = 1100 K at r = 10 mm, convection to 400 K at r = 25 mm, free to expand") {
    const double a = 0.010, b = 0.025, k = 20.0, h = 800.0, Ta = 1100.0, Tinf = 400.0;
    const double C = h * (Tinf - Ta) / (k / b + h * std::log(b / a));
    const verification::ThermoelasticRing ring{a, b, E, alpha, Ta, C};

    struct Result {
      double stress_l2;         // relative L2 error of the stress tensor
      double hoop_inner, hoop_outer;
      double max_reaction;
    };
    // Conduction first, then the nodal temperatures drive the structural model.
    auto run = [&](int n_theta, int n_radial) {
      const Mesh m = meshgen::annulus(a, b, n_theta, n_radial);
      ThermalProblem tp;
      tp.conductivity = k;
      tp.dirichlet.push_back({"inner", [&](const Vec2&) { return Ta; }});
      tp.convection.push_back({"outer", h, Tinf});
      const ThermalSolution ts = solve_thermal(m, tp);

      ElasticProblem ep;
      ep.material = {E, nu, alpha};
      ep.temperature = ts.temperature;
      ep.T_ref = 293.15;
      hold_free(ep, Vec2(a, 0.0), Vec2(-a, 0.0));
      const ElasticSolution es = solve_elastic(m, ep);

      double num = 0.0, den = 0.0;
      std::vector<double> hoop(m.triangle_count());
      for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        const Eigen::Vector3d ex = ring.stress(m.centroid(t));
        num += (es.stress[t] - ex).squaredNorm() * m.area(t);
        den += ex.squaredNorm() * m.area(t);
        const Vec2 c = m.centroid(t).normalized();
        const auto& s = es.stress[t];
        hoop[t] = s[0] * c.y() * c.y() + s[1] * c.x() * c.x() - 2 * s[2] * c.x() * c.y();
      }
      const auto hoop_nodal = nodal_average(m, hoop);
      return Result{std::sqrt(num / den), hoop_nodal[static_cast<std::size_t>(m.nearest_node(Vec2(0.0, a)))],
                    hoop_nodal[static_cast<std::size_t>(m.nearest_node(Vec2(0.0, b)))],
                    es.reactions.cwiseAbs().maxCoeff()};
    };

    WHEN("solved on a 96 x 12 and a 192 x 24 ring mesh") {
      const Result coarse = run(96, 12), fine = run(192, 24);
      const double exact_inner = ring.stress(Vec2(a, 0.0))[1];
      const double exact_outer = ring.stress(Vec2(b, 0.0))[1];

      THEN("the stress field converges to the closed form at first order") {
        CHECK(std::log2(coarse.stress_l2 / fine.stress_l2) == Catch::Approx(1.0).margin(0.1));
        CHECK(fine.stress_l2 < 0.035);
      }
      THEN("the hot inner surface is in hoop compression and the cold outer one in tension") {
        CHECK(fine.hoop_inner < 0.0);
        CHECK(fine.hoop_outer > 0.0);
        // Nodal averaging under-reads boundary peaks; the error still halves with h.
        CHECK(std::abs(fine.hoop_inner / exact_inner - 1.0) < 0.6 * std::abs(coarse.hoop_inner / exact_inner - 1.0));
        CHECK_THAT(fine.hoop_inner, WithinRel(exact_inner, 0.10));
        CHECK_THAT(fine.hoop_outer, WithinRel(exact_outer, 0.05));
      }
      THEN("the supports stay unloaded") { CHECK(fine.max_reaction < 1e-12 * E * alpha * 700.0 * b); }
    }
  }
}
