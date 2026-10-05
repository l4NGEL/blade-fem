#include <bladefem/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace bladefem {

namespace {

constexpr double kPi = 3.14159265358979323846;

double signed_area(const Vec2& a, const Vec2& b, const Vec2& c) {
  return 0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()));
}

void require_cells(int n, const char* what) {
  if (n < 1) throw std::invalid_argument(std::string("meshgen: ") + what + " must be >= 1");
}

// Splits the structured quad (n00, n10, n11, n01) into two triangles. The diagonal
// alternates in a checkerboard so the mesh has no preferred direction.
void add_quad(Mesh& m, int n00, int n10, int n11, int n01, bool flip) {
  if (!flip) {
    m.triangles.push_back({n00, n10, n11});
    m.triangles.push_back({n00, n11, n01});
  } else {
    m.triangles.push_back({n00, n10, n01});
    m.triangles.push_back({n10, n11, n01});
  }
  m.triangle_tags.push_back(0);
  m.triangle_tags.push_back(0);
}

}  // namespace

int Mesh::tag(const std::string& name) const {
  const auto it = names.find(name);
  if (it != names.end()) return it->second;
  std::ostringstream msg;
  msg << "mesh has no boundary named '" << name << "' (known:";
  for (const auto& [n, t] : names) msg << ' ' << n;
  msg << ')';
  throw std::out_of_range(msg.str());
}

std::vector<int> Mesh::boundary_nodes(int t) const {
  std::vector<int> out;
  for (const auto& e : edges) {
    if (e.tag == t) out.insert(out.end(), e.nodes.begin(), e.nodes.end());
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

double Mesh::area(std::size_t t) const {
  const auto& tri = triangles[t];
  return signed_area(nodes[static_cast<std::size_t>(tri[0])], nodes[static_cast<std::size_t>(tri[1])],
                     nodes[static_cast<std::size_t>(tri[2])]);
}

double Mesh::total_area() const {
  double a = 0.0;
  for (std::size_t t = 0; t < triangles.size(); ++t) a += area(t);
  return a;
}

Vec2 Mesh::centroid(std::size_t t) const {
  const auto& tri = triangles[t];
  return (nodes[static_cast<std::size_t>(tri[0])] + nodes[static_cast<std::size_t>(tri[1])] +
          nodes[static_cast<std::size_t>(tri[2])]) / 3.0;
}

int Mesh::nearest_node(const Vec2& p) const {
  if (nodes.empty()) throw std::logic_error("nearest_node on an empty mesh");
  int best = 0;
  double best_d2 = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const double d2 = (nodes[i] - p).squaredNorm();
    if (d2 < best_d2) {
      best_d2 = d2;
      best = static_cast<int>(i);
    }
  }
  return best;
}

double Mesh::max_edge_length() const {
  double h = 0.0;
  for (const auto& tri : triangles) {
    for (int k = 0; k < 3; ++k) {
      const Vec2& a = nodes[static_cast<std::size_t>(tri[static_cast<std::size_t>(k)])];
      const Vec2& b = nodes[static_cast<std::size_t>(tri[static_cast<std::size_t>((k + 1) % 3)])];
      h = std::max(h, (a - b).norm());
    }
  }
  return h;
}

void Mesh::orient_counter_clockwise() {
  for (std::size_t t = 0; t < triangles.size(); ++t) {
    auto& tri = triangles[t];
    const double a = area(t);
    if (a == 0.0) {
      std::ostringstream msg;
      msg << "degenerate (zero-area) triangle " << t;
      throw std::runtime_error(msg.str());
    }
    if (a < 0.0) std::swap(tri[1], tri[2]);
  }
}

namespace meshgen {

Mesh rectangle(double x0, double y0, double x1, double y1, int nx, int ny) {
  require_cells(nx, "nx");
  require_cells(ny, "ny");
  if (!(x1 > x0 && y1 > y0)) throw std::invalid_argument("meshgen::rectangle: empty box");

  Mesh m;
  const auto id = [&](int i, int j) { return j * (nx + 1) + i; };
  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i <= nx; ++i) {
      // Exact end coordinates so that boundary nodes sit exactly on the box.
      const double x = (i == nx) ? x1 : x0 + (x1 - x0) * i / nx;
      const double y = (j == ny) ? y1 : y0 + (y1 - y0) * j / ny;
      m.nodes.emplace_back(x, y);
    }
  }
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      add_quad(m, id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1), (i + j) % 2 == 1);
    }
  }
  m.names = {{"bottom", 1}, {"right", 2}, {"top", 3}, {"left", 4}};
  for (int i = 0; i < nx; ++i) m.edges.push_back({{id(i, 0), id(i + 1, 0)}, 1});
  for (int j = 0; j < ny; ++j) m.edges.push_back({{id(nx, j), id(nx, j + 1)}, 2});
  for (int i = 0; i < nx; ++i) m.edges.push_back({{id(i + 1, ny), id(i, ny)}, 3});
  for (int j = 0; j < ny; ++j) m.edges.push_back({{id(0, j + 1), id(0, j)}, 4});
  return m;
}

Mesh quarter_plate_with_hole(double a, double L, int n_half, int n_radial, double grading) {
  require_cells(n_half, "n_half");
  require_cells(n_radial, "n_radial");
  if (!(a > 0.0 && L > a)) throw std::invalid_argument("meshgen::quarter_plate_with_hole: need 0 < a < L");
  if (!(grading >= 1.0)) throw std::invalid_argument("meshgen::quarter_plate_with_hole: grading must be >= 1");

  // Structured grid: i runs along the hole arc (theta from 0 to 90 deg), j along
  // straight radial rays from the hole to the plate edge.
  const int ni = 2 * n_half;
  Mesh m;
  const auto id = [&](int i, int j) { return j * (ni + 1) + i; };
  for (int j = 0; j <= n_radial; ++j) {
    const double s = std::pow(static_cast<double>(j) / n_radial, grading);
    for (int i = 0; i <= ni; ++i) {
      const double theta = 0.5 * kPi * i / ni;
      Vec2 hole(a * std::cos(theta), a * std::sin(theta));
      Vec2 edge;
      if (i < n_half) {
        edge = Vec2(L, L * std::tan(theta));
      } else if (i == n_half) {
        edge = Vec2(L, L);
      } else {
        edge = Vec2(L / std::tan(theta), L);
      }
      if (i == 0) hole.y() = edge.y() = 0.0;
      if (i == ni) hole.x() = edge.x() = 0.0;
      Vec2 p = hole + s * (edge - hole);
      if (j == n_radial) p = edge;
      m.nodes.push_back(p);
    }
  }
  for (int j = 0; j < n_radial; ++j) {
    for (int i = 0; i < ni; ++i) {
      add_quad(m, id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1), (i + j) % 2 == 1);
    }
  }
  m.orient_counter_clockwise();

  m.names = {{"hole", 1}, {"right", 2}, {"top", 3}, {"symmetry_x", 4}, {"symmetry_y", 5}};
  for (int i = 0; i < ni; ++i) m.edges.push_back({{id(i, 0), id(i + 1, 0)}, 1});
  for (int i = 0; i < ni; ++i) m.edges.push_back({{id(i, n_radial), id(i + 1, n_radial)}, i < n_half ? 2 : 3});
  for (int j = 0; j < n_radial; ++j) {
    m.edges.push_back({{id(ni, j), id(ni, j + 1)}, 4});
    m.edges.push_back({{id(0, j), id(0, j + 1)}, 5});
  }
  return m;
}

Mesh annulus(double r_in, double r_out, int n_theta, int n_radial) {
  if (n_theta < 3) throw std::invalid_argument("meshgen::annulus: n_theta must be >= 3");
  require_cells(n_radial, "n_radial");
  if (!(r_in > 0.0 && r_out > r_in)) throw std::invalid_argument("meshgen::annulus: need 0 < r_in < r_out");

  Mesh m;
  const auto id = [&](int i, int j) { return j * n_theta + (i % n_theta); };
  for (int j = 0; j <= n_radial; ++j) {
    const double r = (j == n_radial) ? r_out : r_in + (r_out - r_in) * j / n_radial;
    for (int i = 0; i < n_theta; ++i) {
      const double theta = 2.0 * kPi * i / n_theta;
      m.nodes.emplace_back(r * std::cos(theta), r * std::sin(theta));
    }
  }
  for (int j = 0; j < n_radial; ++j) {
    for (int i = 0; i < n_theta; ++i) {
      add_quad(m, id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1), (i + j) % 2 == 1);
    }
  }
  m.orient_counter_clockwise();

  m.names = {{"inner", 1}, {"outer", 2}};
  for (int i = 0; i < n_theta; ++i) {
    m.edges.push_back({{id(i + 1, 0), id(i, 0)}, 1});
    m.edges.push_back({{id(i, n_radial), id(i + 1, n_radial)}, 2});
  }
  return m;
}

}  // namespace meshgen

}  // namespace bladefem
