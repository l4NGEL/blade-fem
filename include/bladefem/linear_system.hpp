#pragma once

// Sparse symmetric system K u = f with prescribed (Dirichlet) degrees of freedom.
// Fixed dofs are eliminated, so the reduced matrix stays symmetric positive definite
// and can be factorised with a sparse LDL^T.

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

namespace bladefem {

class LinearSystem {
 public:
  explicit LinearSystem(int dofs);

  int size() const { return n_; }

  void add_matrix(int i, int j, double v) { triplets_.emplace_back(i, j, v); }
  void add_rhs(int i, double v) { f_[i] += v; }

  // Prescribe u_i = value. Repeated calls for the same dof keep the last value.
  void fix(int i, double value);
  bool is_fixed(int i) const { return fixed_[static_cast<std::size_t>(i)]; }

  struct Solution {
    Eigen::VectorXd u;
    // K u - f: the generalised reaction at fixed dofs (zero elsewhere up to round-off).
    Eigen::VectorXd reactions;
  };

  // Throws std::runtime_error if the reduced matrix is singular (e.g. missing constraints).
  Solution solve() const;

 private:
  int n_;
  std::vector<Eigen::Triplet<double>> triplets_;
  Eigen::VectorXd f_;
  std::vector<bool> fixed_;
  Eigen::VectorXd fixed_value_;
};

}  // namespace bladefem
