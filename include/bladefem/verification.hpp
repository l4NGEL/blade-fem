#pragma once

// Error norms and closed-form reference solutions used by the tests and the
// convergence study (apps/bladefem_verify.cpp).

#include <bladefem/mesh.hpp>

#include <Eigen/Core>

#include <functional>
#include <vector>

namespace bladefem::verification {

// ||T_h - T||_L2 for a nodal piecewise-linear field, integrated with a
// 7-point (degree 5) triangle rule.
double l2_error(const Mesh& mesh, const Eigen::VectorXd& nodal,
                const std::function<double(const Vec2&)>& exact);

// |T_h - T|_H1 (gradient seminorm).
double h1_seminorm_error(const Mesh& mesh, const Eigen::VectorXd& nodal,
                         const std::function<Vec2(const Vec2&)>& exact_gradient);

// L2 norm of an exact scalar field, for relative errors.
double l2_norm(const Mesh& mesh, const std::function<double(const Vec2&)>& f);

// Observed order of convergence between two refinement levels.
double observed_rate(double error_coarse, double error_fine, double h_coarse, double h_fine);

// Cantilever of length L, depth D, unit thickness, clamped at x = 0 and loaded
// by a parabolic shear resultant P (downwards) at x = L (Timoshenko & Goodier;
// form used by Augarde & Deeks, Finite Elem. Anal. Des. 44 (2008)).
struct TimoshenkoBeam {
  double L, D, E, nu, P;

  double I() const { return D * D * D / 12.0; }
  Vec2 displacement(const Vec2& p) const;
  Eigen::Vector3d stress(const Vec2& p) const;  // sxx, syy, sxy
  double tip_deflection() const { return displacement(Vec2(L, 0.0)).y(); }
};

// Infinite plate with a circular hole of radius a under remote uniaxial
// tension s0 along x (Kirsch, 1898).
struct KirschPlate {
  double a, s0;
  Eigen::Vector3d stress(const Vec2& p) const;  // sxx, syy, sxy
};

// Free thin ring a <= r <= b (plane stress) carrying the steady conduction
// profile T(r) = T_a + C ln(r / a). Thermal stresses from Timoshenko & Goodier,
// Theory of Elasticity, sec. 151; they depend only on the non-uniform part of T.
struct ThermoelasticRing {
  double a, b, E, alpha, T_a, C;

  double temperature(double r) const;
  // Thin ring (plane stress): sxx, syy, sxy.
  Eigen::Vector3d stress(const Vec2& p) const;
  // Long cylinder with free ends (generalized plane strain, sec. 150):
  // sxx, syy, sxy, szz.
  Eigen::Vector4d stress_long_cylinder(const Vec2& p, double nu) const;
};

}  // namespace bladefem::verification
