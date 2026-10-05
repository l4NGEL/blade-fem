#include <bladefem/elastic2d.hpp>
#include <bladefem/linear_system.hpp>

#include "quadrature.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <string>

namespace bladefem {

namespace {

using Matrix3 = Eigen::Matrix3d;
using Matrix36 = Eigen::Matrix<double, 3, 6>;
using Vector6 = Eigen::Matrix<double, 6, 1>;

Matrix3 constitutive(const ElasticMaterial& mat, PlaneModel model) {
  const double E = mat.E, nu = mat.nu;
  Matrix3 D = Matrix3::Zero();
  if (model == PlaneModel::Stress) {
    const double c = E / (1.0 - nu * nu);
    D << c, c * nu, 0.0, c * nu, c, 0.0, 0.0, 0.0, c * (1.0 - nu) / 2.0;
  } else {
    const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D << c * (1.0 - nu), c * nu, 0.0, c * nu, c * (1.0 - nu), 0.0, 0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0;
  }
  return D;
}

// Strain-displacement matrix of the constant-strain triangle.
Matrix36 strain_matrix(const Vec2& p0, const Vec2& p1, const Vec2& p2, double area) {
  const double b[3] = {p1.y() - p2.y(), p2.y() - p0.y(), p0.y() - p1.y()};
  const double c[3] = {p2.x() - p1.x(), p0.x() - p2.x(), p1.x() - p0.x()};
  Matrix36 B = Matrix36::Zero();
  for (int i = 0; i < 3; ++i) {
    B(0, 2 * i) = b[i];
    B(1, 2 * i + 1) = c[i];
    B(2, 2 * i) = c[i];
    B(2, 2 * i + 1) = b[i];
  }
  return B / (2.0 * area);
}

// Free thermal strain (engineering shear = 0). In plane strain the suppressed
// out-of-plane expansion adds the (1 + nu) factor.
Eigen::Vector3d thermal_strain(const ElasticMaterial& mat, PlaneModel model, double dT) {
  const double e = mat.alpha * dT * (model == PlaneModel::Strain ? 1.0 + mat.nu : 1.0);
  return {e, e, 0.0};
}

void validate(const Mesh& mesh, const ElasticProblem& p) {
  if (!(p.material.E > 0.0)) throw std::invalid_argument("solve_elastic: E must be positive");
  if (!(p.material.nu > -1.0 && p.material.nu < 0.5)) {
    throw std::invalid_argument("solve_elastic: Poisson's ratio must lie in (-1, 0.5)");
  }
  if (p.temperature.size() != 0 && static_cast<std::size_t>(p.temperature.size()) != mesh.node_count()) {
    throw std::invalid_argument("solve_elastic: temperature has " + std::to_string(p.temperature.size()) +
                                " values for " + std::to_string(mesh.node_count()) + " nodes");
  }
  for (const auto& d : p.displacement) {
    if (d.component != 0 && d.component != 1) throw std::invalid_argument("solve_elastic: component must be 0 or 1");
  }
}

}  // namespace

Vec2 ElasticSolution::total_reaction() const {
  Vec2 r = Vec2::Zero();
  for (Eigen::Index i = 0; i + 1 < reactions.size(); i += 2) r += Vec2(reactions[i], reactions[i + 1]);
  return r;
}

double von_mises(const Eigen::Vector3d& s, double szz) {
  const double sx = s[0], sy = s[1], txy = s[2];
  return std::sqrt(0.5 * ((sx - sy) * (sx - sy) + (sy - szz) * (sy - szz) + (szz - sx) * (szz - sx)) +
                   3.0 * txy * txy);
}

std::vector<double> nodal_average(const Mesh& mesh, const std::vector<double>& cell) {
  if (cell.size() != mesh.triangle_count()) throw std::invalid_argument("nodal_average: one value per triangle expected");
  std::vector<double> sum(mesh.node_count(), 0.0), weight(mesh.node_count(), 0.0);
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const double A = mesh.area(t);
    for (int n : mesh.triangles[t]) {
      sum[static_cast<std::size_t>(n)] += A * cell[t];
      weight[static_cast<std::size_t>(n)] += A;
    }
  }
  for (std::size_t i = 0; i < sum.size(); ++i) {
    if (weight[i] > 0.0) sum[i] /= weight[i];
  }
  return sum;
}

ElasticSolution solve_elastic(const Mesh& mesh, const ElasticProblem& p) {
  validate(mesh, p);
  const int n = static_cast<int>(mesh.node_count());
  const bool thermal = p.temperature.size() != 0;
  const Matrix3 D = constitutive(p.material, p.model);
  LinearSystem sys(2 * n);

  const auto element = [&](std::size_t t) {
    const auto& tri = mesh.triangles[t];
    return std::array<int, 6>{2 * tri[0], 2 * tri[0] + 1, 2 * tri[1], 2 * tri[1] + 1, 2 * tri[2], 2 * tri[2] + 1};
  };
  const auto element_dT = [&](std::size_t t) {
    const auto& tri = mesh.triangles[t];
    return (p.temperature[tri[0]] + p.temperature[tri[1]] + p.temperature[tri[2]]) / 3.0 - p.T_ref;
  };
  const auto B_of = [&](std::size_t t) {
    const auto& tri = mesh.triangles[t];
    return strain_matrix(mesh.nodes[static_cast<std::size_t>(tri[0])], mesh.nodes[static_cast<std::size_t>(tri[1])],
                         mesh.nodes[static_cast<std::size_t>(tri[2])], mesh.area(t));
  };

  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const double A = mesh.area(t);
    const Matrix36 B = B_of(t);
    const Eigen::Matrix<double, 6, 6> Ke = A * B.transpose() * D * B;
    const auto dofs = element(t);
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) sys.add_matrix(dofs[static_cast<std::size_t>(i)], dofs[static_cast<std::size_t>(j)], Ke(i, j));
    }
    if (thermal) {
      const Vector6 fe = A * B.transpose() * D * thermal_strain(p.material, p.model, element_dT(t));
      for (int i = 0; i < 6; ++i) sys.add_rhs(dofs[static_cast<std::size_t>(i)], fe[i]);
    }
  }

  for (const auto& tr : p.traction) {
    const int tag = mesh.tag(tr.boundary);
    for (const auto& e : mesh.edges) {
      if (e.tag != tag) continue;
      const Vec2& pa = mesh.nodes[static_cast<std::size_t>(e.nodes[0])];
      const Vec2& pb = mesh.nodes[static_cast<std::size_t>(e.nodes[1])];
      const double l = (pb - pa).norm();
      for (const auto& q : quadrature::kLine3) {
        const Vec2 f = q.w * l * tr.traction((1.0 - q.s) * pa + q.s * pb);
        sys.add_rhs(2 * e.nodes[0], (1.0 - q.s) * f.x());
        sys.add_rhs(2 * e.nodes[0] + 1, (1.0 - q.s) * f.y());
        sys.add_rhs(2 * e.nodes[1], q.s * f.x());
        sys.add_rhs(2 * e.nodes[1] + 1, q.s * f.y());
      }
    }
  }

  for (const auto& d : p.displacement) {
    for (int node : mesh.boundary_nodes(d.boundary)) {
      sys.fix(2 * node + d.component, d.value(mesh.nodes[static_cast<std::size_t>(node)]));
    }
  }
  for (const auto& s : p.supports) {
    const int node = mesh.nearest_node(s.at);
    if (s.fix_x) sys.fix(2 * node, 0.0);
    if (s.fix_y) sys.fix(2 * node + 1, 0.0);
  }

  const LinearSystem::Solution sol = sys.solve();

  ElasticSolution out;
  out.displacement = sol.u;
  out.reactions = sol.reactions;
  out.stress.resize(mesh.triangle_count());
  out.stress_zz.resize(mesh.triangle_count());
  out.von_mises.resize(mesh.triangle_count());
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const auto dofs = element(t);
    Vector6 ue;
    for (int i = 0; i < 6; ++i) ue[i] = sol.u[dofs[static_cast<std::size_t>(i)]];
    const double dT = thermal ? element_dT(t) : 0.0;
    const Eigen::Vector3d sigma = D * (B_of(t) * ue - thermal_strain(p.material, p.model, dT));
    const double szz = p.model == PlaneModel::Strain
                           ? p.material.nu * (sigma[0] + sigma[1]) - p.material.E * p.material.alpha * dT
                           : 0.0;
    out.stress[t] = sigma;
    out.stress_zz[t] = szz;
    out.von_mises[t] = von_mises(sigma, szz);
  }
  return out;
}

}  // namespace bladefem
