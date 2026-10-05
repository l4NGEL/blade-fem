#include <bladefem/linear_system.hpp>
#include <bladefem/thermal2d.hpp>

#include "quadrature.hpp"

#include <stdexcept>

namespace bladefem {

ThermalSolution solve_thermal(const Mesh& mesh, const ThermalProblem& prob) {
  if (!(prob.conductivity > 0.0)) throw std::invalid_argument("solve_thermal: conductivity must be positive");
  const int n = static_cast<int>(mesh.node_count());
  LinearSystem sys(n);

  // Conduction and volumetric source.
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const auto& tri = mesh.triangles[t];
    const Vec2& p0 = mesh.nodes[static_cast<std::size_t>(tri[0])];
    const Vec2& p1 = mesh.nodes[static_cast<std::size_t>(tri[1])];
    const Vec2& p2 = mesh.nodes[static_cast<std::size_t>(tri[2])];
    const double A = mesh.area(t);
    // Constant shape-function gradients of the linear triangle.
    const double b[3] = {p1.y() - p2.y(), p2.y() - p0.y(), p0.y() - p1.y()};
    const double c[3] = {p2.x() - p1.x(), p0.x() - p2.x(), p1.x() - p0.x()};
    const double scale = prob.conductivity / (4.0 * A);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        sys.add_matrix(tri[static_cast<std::size_t>(i)], tri[static_cast<std::size_t>(j)],
                       scale * (b[i] * b[j] + c[i] * c[j]));
      }
    }
    if (prob.source) {
      for (const auto& q : quadrature::kTriangle7) {
        const Vec2 x = q.l0 * p0 + q.l1 * p1 + q.l2 * p2;
        const double wQ = q.w * A * prob.source(x);
        sys.add_rhs(tri[0], wQ * q.l0);
        sys.add_rhs(tri[1], wQ * q.l1);
        sys.add_rhs(tri[2], wQ * q.l2);
      }
    }
  }

  const auto edge_length = [&](const BoundaryEdge& e) {
    return (mesh.nodes[static_cast<std::size_t>(e.nodes[0])] - mesh.nodes[static_cast<std::size_t>(e.nodes[1])]).norm();
  };

  for (const auto& bc : prob.convection) {
    const int tag = mesh.tag(bc.boundary);
    for (const auto& e : mesh.edges) {
      if (e.tag != tag) continue;
      const double l = edge_length(e);
      const int a = e.nodes[0], b = e.nodes[1];
      sys.add_matrix(a, a, bc.h * l / 3.0);
      sys.add_matrix(b, b, bc.h * l / 3.0);
      sys.add_matrix(a, b, bc.h * l / 6.0);
      sys.add_matrix(b, a, bc.h * l / 6.0);
      sys.add_rhs(a, bc.h * bc.T_inf * l / 2.0);
      sys.add_rhs(b, bc.h * bc.T_inf * l / 2.0);
    }
  }
  for (const auto& bc : prob.flux) {
    const int tag = mesh.tag(bc.boundary);
    for (const auto& e : mesh.edges) {
      if (e.tag != tag) continue;
      const double l = edge_length(e);
      sys.add_rhs(e.nodes[0], bc.q_in * l / 2.0);
      sys.add_rhs(e.nodes[1], bc.q_in * l / 2.0);
    }
  }

  // Dirichlet nodes; a node shared by two Dirichlet boundaries is credited to the first.
  std::vector<int> owner(static_cast<std::size_t>(n), -1);
  for (std::size_t k = 0; k < prob.dirichlet.size(); ++k) {
    const auto& bc = prob.dirichlet[k];
    for (int node : mesh.boundary_nodes(bc.boundary)) {
      if (owner[static_cast<std::size_t>(node)] >= 0) continue;
      owner[static_cast<std::size_t>(node)] = static_cast<int>(k);
      sys.fix(node, bc.temperature(mesh.nodes[static_cast<std::size_t>(node)]));
    }
  }

  const LinearSystem::Solution s = sys.solve();

  ThermalSolution out;
  out.temperature = s.u;
  for (std::size_t k = 0; k < prob.dirichlet.size(); ++k) {
    double q = 0.0;
    for (int i = 0; i < n; ++i) {
      if (owner[static_cast<std::size_t>(i)] == static_cast<int>(k)) q += s.reactions[i];
    }
    out.boundary_heat.emplace_back(prob.dirichlet[k].boundary, q);
  }
  for (const auto& bc : prob.convection) {
    const int tag = mesh.tag(bc.boundary);
    double q = 0.0;
    for (const auto& e : mesh.edges) {
      if (e.tag != tag) continue;
      const double Tmean = 0.5 * (s.u[e.nodes[0]] + s.u[e.nodes[1]]);
      q += bc.h * (bc.T_inf - Tmean) * edge_length(e);
    }
    out.boundary_heat.emplace_back(bc.boundary, q);
  }
  for (const auto& bc : prob.flux) {
    const int tag = mesh.tag(bc.boundary);
    double length = 0.0;
    for (const auto& e : mesh.edges) {
      if (e.tag == tag) length += edge_length(e);
    }
    out.boundary_heat.emplace_back(bc.boundary, bc.q_in * length);
  }
  return out;
}

}  // namespace bladefem
