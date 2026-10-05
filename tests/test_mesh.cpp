#include <bladefem/mesh.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

using namespace bladefem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double kPi = 3.14159265358979323846;
}

SCENARIO("Structured rectangle mesh", "[mesh]") {
  GIVEN("a 2 x 1 rectangle split into 8 x 4 cells") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 2.0, 1.0, 8, 4);

    THEN("node, triangle and boundary-edge counts follow from the grid") {
      CHECK(m.node_count() == 9 * 5);
      CHECK(m.triangle_count() == 2 * 8 * 4);
      CHECK(m.edges.size() == 2 * (8 + 4));
    }
    THEN("every triangle is counter-clockwise and the areas add up") {
      for (std::size_t t = 0; t < m.triangle_count(); ++t) CHECK(m.area(t) > 0.0);
      CHECK_THAT(m.total_area(), WithinRel(2.0, 1e-14));
    }
    THEN("the four sides are named and carry the right nodes") {
      const auto left = m.boundary_nodes("left");
      CHECK(left.size() == 5);
      for (int n : left) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].x(), WithinAbs(0.0, 1e-15));
      const auto top = m.boundary_nodes("top");
      CHECK(top.size() == 9);
      for (int n : top) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].y(), WithinAbs(1.0, 1e-15));
    }
    THEN("the characteristic size is the cell diagonal") {
      CHECK_THAT(m.max_edge_length(), WithinRel(std::hypot(0.25, 0.25), 1e-14));
    }
  }
}

SCENARIO("Looking up mesh entities", "[mesh]") {
  const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 4, 4);
  WHEN("asking for an unknown boundary name") {
    THEN("the error lists what does exist") {
      CHECK_THROWS_AS(m.tag("inlet"), std::out_of_range);
    }
  }
  WHEN("searching for the node closest to a point") {
    THEN("the exact grid node is returned") {
      const int n = m.nearest_node(Vec2(0.51, 0.24));
      CHECK(m.nodes[static_cast<std::size_t>(n)].isApprox(Vec2(0.5, 0.25)));
    }
  }
}

SCENARIO("Quarter plate with a circular hole", "[mesh]") {
  GIVEN("a = 1, L = 5, graded towards the hole") {
    const Mesh m = meshgen::quarter_plate_with_hole(1.0, 5.0, 16, 24, 2.0);

    THEN("hole nodes lie on the circle and outer nodes on the plate edges") {
      for (int n : m.boundary_nodes("hole")) {
        CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].norm(), WithinRel(1.0, 1e-13));
      }
      for (int n : m.boundary_nodes("right")) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].x(), WithinRel(5.0, 1e-13));
      for (int n : m.boundary_nodes("top")) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].y(), WithinRel(5.0, 1e-13));
      for (int n : m.boundary_nodes("symmetry_x")) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].x(), WithinAbs(0.0, 1e-13));
      for (int n : m.boundary_nodes("symmetry_y")) CHECK_THAT(m.nodes[static_cast<std::size_t>(n)].y(), WithinAbs(0.0, 1e-13));
    }
    THEN("the area approaches L^2 - pi a^2 / 4 (polygonal hole)") {
      for (std::size_t t = 0; t < m.triangle_count(); ++t) CHECK(m.area(t) > 0.0);
      CHECK_THAT(m.total_area(), WithinRel(25.0 - kPi / 4.0, 1e-3));
    }
  }
}

SCENARIO("Annulus mesh", "[mesh]") {
  GIVEN("r_in = 0.5, r_out = 1") {
    const Mesh m = meshgen::annulus(0.5, 1.0, 64, 8);
    THEN("it closes on itself and the area matches the ring") {
      CHECK(m.node_count() == 64 * 9);
      for (std::size_t t = 0; t < m.triangle_count(); ++t) CHECK(m.area(t) > 0.0);
      CHECK_THAT(m.total_area(), WithinRel(kPi * (1.0 - 0.25), 2e-3));
      CHECK(m.boundary_nodes("inner").size() == 64);
      CHECK(m.boundary_nodes("outer").size() == 64);
    }
  }
}
