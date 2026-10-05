#include <bladefem/mesh.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>
#include <stdexcept>

using namespace bladefem;
using Catch::Matchers::WithinRel;

namespace {

// Unit square, two triangles (the second written clockwise on purpose),
// a "wall" line on x = 0, a "far" line on x = 1 and a vertex element that the
// reader must skip.
constexpr const char* kUnitSquare = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
3
1 1 "wall"
1 2 "far"
2 10 "solid"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
6
1 15 2 99 1 1
2 1 2 1 4 4 1
3 1 2 2 2 2 3
4 2 2 10 1 1 2 3
5 2 2 10 1 1 4 3
6 1 2 99 9 1 2
$EndElements
)";

}  // namespace

SCENARIO("Reading a Gmsh 2.2 ASCII mesh", "[gmsh]") {
  GIVEN("a unit square with named boundaries") {
    std::istringstream in(kUnitSquare);

    WHEN("it is read") {
      const Mesh m = read_gmsh(in);

      THEN("nodes and triangles are converted to zero-based indices") {
        CHECK(m.node_count() == 4);
        CHECK(m.triangle_count() == 2);
        CHECK(m.triangle_tags[0] == 10);
      }
      THEN("clockwise triangles are re-oriented counter-clockwise") {
        for (std::size_t t = 0; t < m.triangle_count(); ++t) CHECK(m.area(t) > 0.0);
        CHECK_THAT(m.total_area(), WithinRel(1.0, 1e-15));
      }
      THEN("physical names map to the boundary edges") {
        CHECK(m.tag("wall") == 1);
        CHECK(m.boundary_nodes("wall") == std::vector<int>{0, 3});
        CHECK(m.boundary_nodes("far") == std::vector<int>{1, 2});
      }
      THEN("unnamed physical groups are still addressable by tag") {
        CHECK(m.boundary_nodes(99) == std::vector<int>{0, 1});
      }
    }
  }

  GIVEN("a mesh in the MSH 4 format") {
    std::istringstream in("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
    THEN("the reader refuses it and says how to convert") {
      CHECK_THROWS_AS(read_gmsh(in), std::runtime_error);
    }
  }

  GIVEN("a path that does not exist") {
    THEN("reading fails with a clear error") {
      CHECK_THROWS_AS(read_gmsh(std::string("does/not/exist.msh")), std::runtime_error);
    }
  }
}
