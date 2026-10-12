#include <bladefem/linear_system.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bladefem {

LinearSystem::LinearSystem(int dofs)
    : n_(dofs),
      f_(Eigen::VectorXd::Zero(dofs)),
      fixed_(static_cast<std::size_t>(dofs), false),
      fixed_value_(Eigen::VectorXd::Zero(dofs)) {
  if (dofs < 1) throw std::invalid_argument("LinearSystem: need at least one dof");
}

void LinearSystem::fix(int i, double value) {
  fixed_[static_cast<std::size_t>(i)] = true;
  fixed_value_[i] = value;
}

LinearSystem::Solution LinearSystem::solve() const {
  Eigen::SparseMatrix<double> K(n_, n_);
  K.setFromTriplets(triplets_.begin(), triplets_.end());

  // Number the free dofs.
  std::vector<int> free_index(static_cast<std::size_t>(n_), -1);
  int nf = 0;
  for (int i = 0; i < n_; ++i) {
    if (!fixed_[static_cast<std::size_t>(i)]) free_index[static_cast<std::size_t>(i)] = nf++;
  }

  Solution out;
  out.u = fixed_value_;

  if (nf > 0) {
    std::vector<Eigen::Triplet<double>> reduced;
    reduced.reserve(static_cast<std::size_t>(K.nonZeros()));
    Eigen::VectorXd rhs(nf);
    for (int i = 0; i < n_; ++i) {
      const int fi = free_index[static_cast<std::size_t>(i)];
      if (fi >= 0) rhs[fi] = f_[i];
    }
    for (int col = 0; col < K.outerSize(); ++col) {
      const int fc = free_index[static_cast<std::size_t>(col)];
      for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
        const int fr = free_index[static_cast<std::size_t>(it.row())];
        if (fr < 0) continue;
        if (fc >= 0) {
          reduced.emplace_back(fr, fc, it.value());
        } else {
          rhs[fr] -= it.value() * fixed_value_[col];  // lift the prescribed value
        }
      }
    }
    Eigen::SparseMatrix<double> Kff(nf, nf);
    Kff.setFromTriplets(reduced.begin(), reduced.end());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(Kff);
    if (ldlt.info() != Eigen::Success) throw std::runtime_error("LinearSystem: factorisation failed");

    // A rigid-body mode or a free temperature level shows up as a (near) zero pivot.
    const Eigen::VectorXd d = ldlt.vectorD();
    const double dmax = d.cwiseAbs().maxCoeff();
    const double dmin = d.minCoeff();
    if (!(dmin > 1e-12 * dmax)) {
      throw std::runtime_error(
          "LinearSystem: matrix is singular or indefinite - the boundary conditions do not "
          "remove all rigid-body / constant modes");
    }
    const Eigen::VectorXd uf = ldlt.solve(rhs);
    for (int i = 0; i < n_; ++i) {
      const int fi = free_index[static_cast<std::size_t>(i)];
      if (fi >= 0) out.u[i] = uf[fi];
    }
  }

  out.reactions = K * out.u - f_;
  for (int i = 0; i < n_; ++i) {
    if (!fixed_[static_cast<std::size_t>(i)]) out.reactions[i] = 0.0;
  }
  return out;
}

}  // namespace bladefem
