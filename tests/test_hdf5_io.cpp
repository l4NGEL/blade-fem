#include <bladefem/hdf5_io.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

using namespace bladefem;
namespace fs = std::filesystem;

SCENARIO("Writing results to HDF5 with an XDMF descriptor", "[hdf5]") {
  GIVEN("a small mesh with a nodal scalar, a nodal vector and a cell scalar") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 2, 1);  // 6 nodes, 4 triangles
    std::vector<double> T(m.node_count()), u(2 * m.node_count()), vm(m.triangle_count());
    for (std::size_t i = 0; i < m.node_count(); ++i) {
      T[i] = 300.0 + static_cast<double>(i);
      u[2 * i] = 0.1 * static_cast<double>(i);
      u[2 * i + 1] = -0.2 * static_cast<double>(i);
    }
    for (std::size_t t = 0; t < m.triangle_count(); ++t) vm[t] = 1.0e6 * static_cast<double>(t + 1);

    const fs::path dir = fs::temp_directory_path() / "bladefem_test_hdf5" / "nested";
    fs::remove_all(dir.parent_path());
    const std::string base = (dir / "case").string();

    WHEN("they are written") {
      write_results(base, m,
                    {{"temperature", Field::Location::Node, 1, T},
                     {"displacement", Field::Location::Node, 2, u},
                     {"von_mises", Field::Location::Cell, 1, vm}});

      THEN("the mesh and fields read back unchanged") {
        const Dataset nodes = read_dataset(base + ".h5", "/mesh/nodes");
        CHECK(nodes.shape == std::vector<std::size_t>{6, 2});
        CHECK(nodes.values[2 * 5] == 1.0);  // node 5 = (1, 1)
        const Dataset tris = read_dataset(base + ".h5", "/mesh/triangles");
        CHECK(tris.shape == std::vector<std::size_t>{4, 3});
        const Dataset temp = read_dataset(base + ".h5", "/fields/temperature");
        CHECK(temp.values == T);
        const Dataset disp = read_dataset(base + ".h5", "/fields/displacement");
        CHECK(disp.shape == std::vector<std::size_t>{6, 3});
        CHECK(disp.values[3 * 4 + 1] == u[2 * 4 + 1]);
        CHECK(disp.values[3 * 4 + 2] == 0.0);
        CHECK(read_dataset(base + ".h5", "/fields/von_mises").values == vm);
      }
      THEN("the XDMF file points at those datasets") {
        std::ifstream x(base + ".xdmf");
        REQUIRE(x.good());
        std::stringstream ss;
        ss << x.rdbuf();
        const std::string xml = ss.str();
        using Catch::Matchers::ContainsSubstring;
        CHECK_THAT(xml, ContainsSubstring("TopologyType=\"Triangle\""));
        CHECK_THAT(xml, ContainsSubstring("case.h5:/mesh/nodes"));
        CHECK_THAT(xml, ContainsSubstring("Name=\"displacement\" AttributeType=\"Vector\" Center=\"Node\""));
        CHECK_THAT(xml, ContainsSubstring("Name=\"von_mises\" AttributeType=\"Scalar\" Center=\"Cell\""));
      }
    }
  }

  GIVEN("a field whose size does not match the mesh") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, 1, 1);
    THEN("writing is refused before anything is created") {
      const std::string base = (fs::temp_directory_path() / "bladefem_bad").string();
      CHECK_THROWS_AS(write_results(base, m, {{"T", Field::Location::Node, 1, {1.0, 2.0}}}), std::invalid_argument);
    }
  }

  GIVEN("a dataset that does not exist") {
    THEN("reading it fails") {
      CHECK_THROWS(read_dataset("no_such_file.h5", "/mesh/nodes"));
    }
  }
}
