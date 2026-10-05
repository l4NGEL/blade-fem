#pragma once

// Quadrature rules shared by the element routines (internal header).

#include <array>

namespace bladefem::quadrature {

// Point on the reference triangle in barycentric coordinates, weight normalised
// so that the weights sum to one (multiply by the element area).
struct TrianglePoint {
  double l0, l1, l2, w;
};

// Dunavant degree-5 rule (7 points).
inline constexpr std::array<TrianglePoint, 7> kTriangle7{{
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225},
    {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506},
    {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506},
    {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506},
    {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
    {0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827},
    {0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827},
}};

// Degree-2 rule (3 interior points).
inline constexpr std::array<TrianglePoint, 3> kTriangle3{{
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0},
    {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0},
}};

// Gauss-Legendre on [0, 1] (3 points, exact to degree 5), weights sum to one.
struct LinePoint {
  double s, w;
};
inline constexpr std::array<LinePoint, 3> kLine3{{
    {0.11270166537925831, 5.0 / 18.0},
    {0.5, 8.0 / 18.0},
    {0.88729833462074169, 5.0 / 18.0},
}};

}  // namespace bladefem::quadrature
