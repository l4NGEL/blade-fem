#include <bladefem/elastic2d.hpp>
#include <bladefem/linear_system.hpp>

#include "quadrature.hpp"

#include <Eigen/Dense>

#include <array>
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
  return std::sqrt(0.5 * ((sx - sy) * (sx - sy) + (sy - szz) * (sy - szz) + (szz - sx) * (szz - sx)) + 3.0 * txy * txy);
}

std::vector<double> nodal_average(const Mesh& mesh, const std::vector<double>& cell) {
  if (cell.size() != mesh.triangle_count())
    throw std::invalid_argument("nodal_average: one value per triangle expected");
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
  const bool gps = p.model == PlaneModel::GeneralizedStrain;
  const double E = p.material.E, nu = p.material.nu, alpha = p.material.alpha;

  // Generalized plane strain appends three section-wide dofs (e0, kx, ky) after
  // the 2N nodal displacements; kx, ky are measured from the area centroid.
  const int g = 2 * n;
  LinearSystem sys(gps ? g + 3 : g);
  Vec2 origin = Vec2::Zero();
  if (gps) {
    for (std::size_t t = 0; t < mesh.triangle_count(); ++t) origin += mesh.area(t) * mesh.centroid(t);
    origin /= mesh.total_area();
  }

  // 3D isotropic stiffness in (xx, yy, xy, zz) order for generalized plane strain.
  Eigen::Matrix4d D4 = Eigen::Matrix4d::Zero();
  if (gps) {
    const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D4 << c * (1 - nu), c * nu, 0, c * nu, c * nu, c * (1 - nu), 0, c * nu, 0, 0, c * (1 - 2 * nu) / 2, 0, c * nu,
        c * nu, 0, c * (1 - nu);
  }
  const Matrix3 D = gps ? Matrix3::Zero() : constitutive(p.material, p.model);

  const auto element_dofs = [&](std::size_t t) {
    const auto& tri = mesh.triangles[t];
    std::vector<int> d{2 * tri[0], 2 * tri[0] + 1, 2 * tri[1], 2 * tri[1] + 1, 2 * tri[2], 2 * tri[2] + 1};
    if (gps) d.insert(d.end(), {g, g + 1, g + 2});
    return d;
  };
  const auto element_dT = [&](std::size_t t) {
    if (!thermal) return 0.0;
    const auto& tri = mesh.triangles[t];
    return (p.temperature[tri[0]] + p.temperature[tri[1]] + p.temperature[tri[2]]) / 3.0 - p.T_ref;
  };
  const auto B_of = [&](std::size_t t) {
    const auto& tri = mesh.triangles[t];
    return strain_matrix(mesh.nodes[static_cast<std::size_t>(tri[0])], mesh.nodes[static_cast<std::size_t>(tri[1])],
                         mesh.nodes[static_cast<std::size_t>(tri[2])], mesh.area(t));
  };
  // Strain operator of the generalized plane strain element at point x.
  const auto B4_of = [&](const Matrix36& B, const Vec2& x) {
    Eigen::Matrix<double, 4, 9> B4 = Eigen::Matrix<double, 4, 9>::Zero();
    B4.topLeftCorner<3, 6>() = B;
    B4(3, 6) = 1.0;
    B4(3, 7) = x.x() - origin.x();
    B4(3, 8) = x.y() - origin.y();
    return B4;
  };

  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const double A = mesh.area(t);
    const Matrix36 B = B_of(t);
    const double dT = element_dT(t);
    Eigen::MatrixXd Ke;
    Eigen::VectorXd fe;
    if (!gps) {
      Ke = A * B.transpose() * D * B;
      fe = A * B.transpose() * D * thermal_strain(p.material, p.model, dT);
    } else {
      // eps_zz varies linearly over the element, so integrate with a degree-2 rule.
      Ke = Eigen::MatrixXd::Zero(9, 9);
      fe = Eigen::VectorXd::Zero(9);
      const Eigen::Vector4d eps0(alpha * dT, alpha * dT, 0.0, alpha * dT);
      const auto& tri = mesh.triangles[t];
      for (const auto& q : quadrature::kTriangle3) {
        const Vec2 x = q.l0 * mesh.nodes[static_cast<std::size_t>(tri[0])] +
                       q.l1 * mesh.nodes[static_cast<std::size_t>(tri[1])] +
                       q.l2 * mesh.nodes[static_cast<std::size_t>(tri[2])];
        const auto B4 = B4_of(B, x);
        Ke += q.w * A * B4.transpose() * D4 * B4;
        fe += q.w * A * B4.transpose() * D4 * eps0;
      }
    }
    const auto dofs = element_dofs(t);
    for (std::size_t i = 0; i < dofs.size(); ++i) {
      for (std::size_t j = 0; j < dofs.size(); ++j) {
        sys.add_matrix(dofs[i], dofs[j], Ke(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
      }
      if (thermal) sys.add_rhs(dofs[i], fe[static_cast<Eigen::Index>(i)]);
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
  out.displacement = sol.u.head(g);
  out.reactions = sol.reactions.head(g);
  if (gps) {
    out.axial_strain = sol.u.tail(3);
    out.axial_origin = origin;
  }
  out.stress.resize(mesh.triangle_count());
  out.stress_zz.resize(mesh.triangle_count());
  out.von_mises.resize(mesh.triangle_count());
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const auto dofs = element_dofs(t);
    Eigen::VectorXd ue(static_cast<Eigen::Index>(dofs.size()));
    for (std::size_t i = 0; i < dofs.size(); ++i) ue[static_cast<Eigen::Index>(i)] = sol.u[dofs[i]];
    const double dT = element_dT(t);
    Eigen::Vector3d sigma;
    double szz = 0.0;
    if (gps) {
      const Eigen::Vector4d eps0(alpha * dT, alpha * dT, 0.0, alpha * dT);
      const Eigen::Vector4d s4 = D4 * (B4_of(B_of(t), mesh.centroid(t)) * ue - eps0);
      sigma = s4.head<3>();
      szz = s4[3];
    } else {
      sigma = D * (B_of(t) * ue - thermal_strain(p.material, p.model, dT));
      if (p.model == PlaneModel::Strain) szz = nu * (sigma[0] + sigma[1]) - E * alpha * dT;
    }
    out.stress[t] = sigma;
    out.stress_zz[t] = szz;
    out.von_mises[t] = von_mises(sigma, szz);
  }
  return out;
}

}  // namespace bladefem
