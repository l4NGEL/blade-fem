#pragma once

// Two-dimensional linear-triangle mesh with tagged boundary edges.

#include <Eigen/Core>

#include <array>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace bladefem {

using Vec2 = Eigen::Vector2d;

struct BoundaryEdge {
  std::array<int, 2> nodes;
  int tag;
};

struct Mesh {
  std::vector<Vec2> nodes;
  std::vector<std::array<int, 3>> triangles;  // counter-clockwise node order
  std::vector<int> triangle_tags;
  std::vector<BoundaryEdge> edges;
  std::map<std::string, int> names;           // physical/boundary name -> tag

  std::size_t node_count() const { return nodes.size(); }
  std::size_t triangle_count() const { return triangles.size(); }

  // Tag registered under `name`; throws std::out_of_range with the known names otherwise.
  int tag(const std::string& name) const;
  bool has(const std::string& name) const { return names.count(name) != 0; }

  // Sorted, de-duplicated nodes lying on the edges carrying `tag`.
  std::vector<int> boundary_nodes(int tag) const;
  std::vector<int> boundary_nodes(const std::string& name) const { return boundary_nodes(tag(name)); }

  double area(std::size_t triangle) const;
  double total_area() const;
  Vec2 centroid(std::size_t triangle) const;
  int nearest_node(const Vec2& p) const;
  // Longest edge over all triangles: the usual "h" in convergence studies.
  double max_edge_length() const;

  // Reorders every triangle counter-clockwise; throws on degenerate (zero-area) ones.
  void orient_counter_clockwise();
};

namespace meshgen {

// [x0,x1] x [y0,y1] split into nx*ny cells, each cut into two triangles with
// alternating diagonals. Boundary names: bottom, right, top, left.
Mesh rectangle(double x0, double y0, double x1, double y1, int nx, int ny);

// Quarter of a square plate [0,L]^2 with a circular hole of radius a at the origin.
// `n_half` cells along each half of the hole arc, `n_radial` cells from hole to edge,
// radial spacing graded towards the hole by `grading` >= 1 (1 = uniform).
// Boundary names: hole, right (x = L), top (y = L), symmetry_x (x = 0), symmetry_y (y = 0).
Mesh quarter_plate_with_hole(double a, double L, int n_half, int n_radial, double grading = 1.0);

// Full annulus r_in <= r <= r_out. Boundary names: inner, outer.
Mesh annulus(double r_in, double r_out, int n_theta, int n_radial);

}  // namespace meshgen

// Gmsh MSH 2.2 ASCII reader (`gmsh -format msh22`). 3-node triangles become the
// mesh, 2-node lines become tagged boundary edges; physical names are kept.
Mesh read_gmsh(const std::string& path);
Mesh read_gmsh(std::istream& in);

}  // namespace bladefem
