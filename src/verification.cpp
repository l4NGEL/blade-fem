#include <bladefem/verification.hpp>

#include "quadrature.hpp"

#include <cmath>
#include <stdexcept>

namespace bladefem::verification {

namespace {

template <class Integrand>
double integrate(const Mesh& mesh, Integrand f) {
  double sum = 0.0;
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const auto& tri = mesh.triangles[t];
    const Vec2& p0 = mesh.nodes[static_cast<std::size_t>(tri[0])];
    const Vec2& p1 = mesh.nodes[static_cast<std::size_t>(tri[1])];
    const Vec2& p2 = mesh.nodes[static_cast<std::size_t>(tri[2])];
    const double A = mesh.area(t);
    for (const auto& q : quadrature::kTriangle7) {
      sum += q.w * A * f(t, q.l0 * p0 + q.l1 * p1 + q.l2 * p2, q);
    }
  }
  return sum;
}

}  // namespace

double l2_error(const Mesh& mesh, const Eigen::VectorXd& u, const std::function<double(const Vec2&)>& exact) {
  const double e2 = integrate(mesh, [&](std::size_t t, const Vec2& x, const quadrature::TrianglePoint& q) {
    const auto& tri = mesh.triangles[t];
    const double uh = q.l0 * u[tri[0]] + q.l1 * u[tri[1]] + q.l2 * u[tri[2]];
    const double e = uh - exact(x);
    return e * e;
  });
  return std::sqrt(e2);
}

double h1_seminorm_error(const Mesh& mesh, const Eigen::VectorXd& u,
                         const std::function<Vec2(const Vec2&)>& exact_gradient) {
  const double e2 = integrate(mesh, [&](std::size_t t, const Vec2& x, const quadrature::TrianglePoint&) {
    const auto& tri = mesh.triangles[t];
    const Vec2& p0 = mesh.nodes[static_cast<std::size_t>(tri[0])];
    const Vec2& p1 = mesh.nodes[static_cast<std::size_t>(tri[1])];
    const Vec2& p2 = mesh.nodes[static_cast<std::size_t>(tri[2])];
    const double twoA = 2.0 * mesh.area(t);
    const Vec2 grad((u[tri[0]] * (p1.y() - p2.y()) + u[tri[1]] * (p2.y() - p0.y()) + u[tri[2]] * (p0.y() - p1.y())) / twoA,
                    (u[tri[0]] * (p2.x() - p1.x()) + u[tri[1]] * (p0.x() - p2.x()) + u[tri[2]] * (p1.x() - p0.x())) / twoA);
    return (grad - exact_gradient(x)).squaredNorm();
  });
  return std::sqrt(e2);
}

double l2_norm(const Mesh& mesh, const std::function<double(const Vec2&)>& f) {
  return std::sqrt(integrate(mesh, [&](std::size_t, const Vec2& x, const quadrature::TrianglePoint&) {
    const double v = f(x);
    return v * v;
  }));
}

double observed_rate(double e_coarse, double e_fine, double h_coarse, double h_fine) {
  if (!(e_coarse > 0.0 && e_fine > 0.0 && h_coarse > h_fine && h_fine > 0.0)) {
    throw std::invalid_argument("observed_rate: need positive errors and h_coarse > h_fine > 0");
  }
  return std::log(e_coarse / e_fine) / std::log(h_coarse / h_fine);
}

Vec2 TimoshenkoBeam::displacement(const Vec2& p) const {
  const double x = p.x(), y = p.y();
  const double c = P / (6.0 * E * I());
  const double ux = c * y * ((6.0 * L - 3.0 * x) * x + (2.0 + nu) * (y * y - D * D / 4.0));
  const double uy = -c * (3.0 * nu * y * y * (L - x) + (4.0 + 5.0 * nu) * D * D * x / 4.0 + (3.0 * L - x) * x * x);
  return {ux, uy};
}

Eigen::Vector3d TimoshenkoBeam::stress(const Vec2& p) const {
  const double x = p.x(), y = p.y();
  return {P * (L - x) * y / I(), 0.0, -P / (2.0 * I()) * (D * D / 4.0 - y * y)};
}

Eigen::Vector3d KirschPlate::stress(const Vec2& p) const {
  const double r2 = p.squaredNorm();
  const double th = std::atan2(p.y(), p.x());
  const double a2 = a * a / r2;
  const double a4 = a2 * a2;
  const double c2 = std::cos(2 * th), c4 = std::cos(4 * th), s2 = std::sin(2 * th), s4 = std::sin(4 * th);
  return {s0 * (1.0 - a2 * (1.5 * c2 + c4) + 1.5 * a4 * c4),
          s0 * (-a2 * (0.5 * c2 - c4) - 1.5 * a4 * c4),
          s0 * (-a2 * (0.5 * s2 + s4) + 1.5 * a4 * s4)};
}

}  // namespace bladefem::verification
