#include <bladefem/fem1d.hpp>

#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace bladefem::fem1d {

EndCondition EndCondition::temperature(double T) {
  return {Kind::Temperature, T, 0.0};
}
EndCondition EndCondition::heat_flux(double q_in) {
  return {Kind::HeatFlux, q_in, 0.0};
}
EndCondition EndCondition::convection(double h, double T_inf) {
  return {Kind::Convection, T_inf, h};
}

double BarSolution::at(double xq) const {
  if (x.empty()) throw std::logic_error("BarSolution::at on empty solution");
  if (xq <= x.front()) return temperature.front();
  if (xq >= x.back()) return temperature.back();
  const auto it = std::upper_bound(x.begin(), x.end(), xq);
  const auto i = static_cast<std::size_t>(it - x.begin()) - 1;
  const double t = (xq - x[i]) / (x[i + 1] - x[i]);
  return (1.0 - t) * temperature[i] + t * temperature[i + 1];
}

namespace {

// 3-point Gauss-Legendre rule on [-1, 1]; exact for polynomials up to degree 5.
constexpr std::array<double, 3> kGaussPoints{-0.7745966692414834, 0.0, 0.7745966692414834};
constexpr std::array<double, 3> kGaussWeights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

void validate(const BarProblem& p) {
  if (p.elements < 1) throw std::invalid_argument("fem1d: need at least one element");
  if (!(p.length > 0.0)) throw std::invalid_argument("fem1d: length must be positive");
  if (!(p.conductivity > 0.0)) throw std::invalid_argument("fem1d: conductivity must be positive");
  const auto anchors = [](const EndCondition& c) {
    return c.kind == EndCondition::Kind::Temperature || (c.kind == EndCondition::Kind::Convection && c.h > 0.0);
  };
  if (!anchors(p.left) && !anchors(p.right)) {
    throw std::invalid_argument(
        "fem1d: at least one end needs a temperature or convection condition, "
        "otherwise the temperature level is undetermined");
  }
}

}  // namespace

BarSolution solve(const BarProblem& p) {
  validate(p);

  const int ne = p.elements;
  const int nn = ne + 1;
  const double he = p.length / ne;

  BarSolution out;
  out.x.resize(static_cast<std::size_t>(nn));
  for (int i = 0; i < nn; ++i) out.x[static_cast<std::size_t>(i)] = i * he;

  Eigen::SparseMatrix<double> K(nn, nn);
  Eigen::VectorXd f = Eigen::VectorXd::Zero(nn);
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(static_cast<std::size_t>(4 * ne + 2));

  const double ke = p.conductivity / he;
  for (int e = 0; e < ne; ++e) {
    const int a = e, b = e + 1;
    trip.emplace_back(a, a, ke);
    trip.emplace_back(a, b, -ke);
    trip.emplace_back(b, a, -ke);
    trip.emplace_back(b, b, ke);

    if (p.source) {
      const double x0 = out.x[static_cast<std::size_t>(a)];
      for (std::size_t g = 0; g < kGaussPoints.size(); ++g) {
        const double xi = kGaussPoints[g];
        const double w = kGaussWeights[g] * he / 2.0;
        const double Q = p.source(x0 + (xi + 1.0) * he / 2.0);
        f[a] += w * Q * (1.0 - xi) / 2.0;
        f[b] += w * Q * (1.0 + xi) / 2.0;
      }
    }
  }

  // Natural (flux, convection) conditions enter the weak form directly.
  const auto apply_natural = [&](const EndCondition& c, int node) {
    switch (c.kind) {
      case EndCondition::Kind::HeatFlux:
        f[node] += c.value;
        break;
      case EndCondition::Kind::Convection:
        trip.emplace_back(node, node, c.h);
        f[node] += c.h * c.value;
        break;
      case EndCondition::Kind::Temperature:
        break;
    }
  };
  apply_natural(p.left, 0);
  apply_natural(p.right, nn - 1);
  K.setFromTriplets(trip.begin(), trip.end());

  // Essential conditions: replace the row by the identity and lift the column
  // into the right-hand side, which keeps the matrix symmetric.
  Eigen::VectorXd fixed_value = Eigen::VectorXd::Zero(nn);
  std::vector<bool> fixed(static_cast<std::size_t>(nn), false);
  if (p.left.kind == EndCondition::Kind::Temperature) {
    fixed[0] = true;
    fixed_value[0] = p.left.value;
  }
  if (p.right.kind == EndCondition::Kind::Temperature) {
    fixed[static_cast<std::size_t>(nn - 1)] = true;
    fixed_value[nn - 1] = p.right.value;
  }
  f -= K * fixed_value;
  for (int k = 0; k < K.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
      if (fixed[static_cast<std::size_t>(it.row())] || fixed[static_cast<std::size_t>(it.col())]) {
        it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
      }
    }
  }
  for (int i = 0; i < nn; ++i) {
    if (fixed[static_cast<std::size_t>(i)]) f[i] = fixed_value[i];
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(K);
  if (solver.info() != Eigen::Success) throw std::runtime_error("fem1d: factorisation failed");
  const Eigen::VectorXd T = solver.solve(f);

  out.temperature.assign(T.data(), T.data() + T.size());
  return out;
}

}  // namespace bladefem::fem1d
