#include <bladefem/case_config.hpp>
#include <bladefem/case_runner.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <sstream>
#include <stdexcept>

using namespace bladefem;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace fs = std::filesystem;

namespace {

constexpr const char* kCase = R"({
  "name": "hot plate",
  "mesh": "meshes/plate.msh",
  "output": "out/plate",
  "thermal": {
    "conductivity": 20.0,
    "temperature": [ { "boundary": "left", "T": 900.0 } ],
    "convection":  [ { "boundary": "right", "h": 500.0, "T_inf": 300.0 } ],
    "flux":        [ { "boundary": "top", "q_in": 0.0 } ]
  },
  "structural": {
    "E": 2.0e11, "nu": 0.3, "alpha": 1.3e-5,
    "model": "generalized_plane_strain",
    "T_ref": 300.0,
    "supports": [ { "at": [0.0, 0.0] }, { "at": [0.1, 0.0], "x": false } ],
    "displacement": [ { "boundary": "bottom", "uy": 0.0 } ],
    "traction": [ { "boundary": "top", "t": [0.0, -1.0e6] } ]
  }
})";

CaseConfig parse(const std::string& json) {
  std::istringstream in(json);
  return parse_case(in, "/data/cases");
}

std::string replace(std::string s, const std::string& from, const std::string& to) {
  const auto pos = s.find(from);
  REQUIRE(pos != std::string::npos);
  return s.replace(pos, from.size(), to);
}

}  // namespace

SCENARIO("Reading an analysis case from JSON", "[case]") {
  GIVEN("a complete thermal + structural case") {
    const CaseConfig c = parse(kCase);

    THEN("paths are resolved against the case file's directory") {
      CHECK(fs::path(c.mesh) == fs::path("/data/cases/meshes/plate.msh"));
      CHECK(fs::path(c.output) == fs::path("/data/cases/out/plate"));
    }
    THEN("every thermal boundary condition is read") {
      REQUIRE(c.thermal);
      CHECK(c.thermal->conductivity == 20.0);
      REQUIRE(c.thermal->convection.size() == 1);
      CHECK(c.thermal->convection[0].h == 500.0);
      CHECK(c.thermal->temperature[0].T == 900.0);
      CHECK(c.thermal->flux[0].boundary == "top");
    }
    THEN("the structural model, supports and loads are read") {
      REQUIRE(c.structural);
      CHECK(c.structural->model == PlaneModel::GeneralizedStrain);
      CHECK(c.structural->material.alpha == 1.3e-5);
      REQUIRE(c.structural->supports.size() == 2);
      CHECK(c.structural->supports[1].x == false);
      CHECK(c.structural->supports[1].y == true);
      CHECK(c.structural->supports[1].at.isApprox(Vec2(0.1, 0.0)));
      CHECK_FALSE(c.structural->displacement[0].ux.has_value());
      CHECK(*c.structural->displacement[0].uy == 0.0);
      CHECK(c.structural->traction[0].traction.y() == -1.0e6);
    }
  }

  GIVEN("a case without a structural section") {
    const std::string json = R"({"mesh": "m.msh", "thermal": {"conductivity": 1.0}})";
    THEN("only the thermal analysis is configured") {
      const CaseConfig c = parse(json);
      CHECK(c.thermal);
      CHECK_FALSE(c.structural);
      CHECK(c.output.empty());
    }
  }

  GIVEN("mistakes in the file") {
    THEN("a missing material property is named in the error") {
      CHECK_THROWS_WITH(parse(replace(kCase, R"("conductivity": 20.0,)", "")),
                        ContainsSubstring("thermal.conductivity"));
    }
    THEN("a misspelt model lists the valid choices") {
      CHECK_THROWS_WITH(parse(replace(kCase, "generalized_plane_strain", "plane_strian")),
                        ContainsSubstring("plane_stress"));
    }
    THEN("a non-numeric value is reported") {
      CHECK_THROWS_WITH(parse(replace(kCase, R"("h": 500.0)", R"("h": "high")")),
                        ContainsSubstring("thermal.convection[0].h"));
    }
    THEN("a case with neither analysis is rejected") {
      CHECK_THROWS_AS(parse(R"({"mesh": "m.msh"})"), std::runtime_error);
    }
    THEN("malformed JSON is reported as such") {
      CHECK_THROWS_WITH(parse("{ \"mesh\": "), ContainsSubstring("JSON"));
    }
  }
}

SCENARIO("Running a coupled case end to end", "[case][integration]") {
  GIVEN("a plate held at 900 K on the left and cooled on the right, free to expand") {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 0.1, 0.05, 20, 10);
    CaseConfig c = parse(kCase);
    c.structural->displacement.clear();
    c.structural->traction.clear();
    c.output = (fs::temp_directory_path() / "bladefem_case_test" / "plate").string();
    fs::remove_all(fs::path(c.output).parent_path());

    WHEN("it is run") {
      const CaseReport r = run_case(c, m);

      THEN("temperatures lie between the coolant and the hot wall, and heat is conserved") {
        CHECK(r.thermal);
        CHECK_THAT(r.T_max, WithinRel(900.0, 1e-12));
        CHECK(r.T_min > 300.0);
        double net = 0.0;
        for (const auto& [name, q] : r.boundary_heat) net += q;
        CHECK_THAT(net, WithinAbs(0.0, 1e-9 * std::abs(r.boundary_heat.front().second)));
      }
      THEN("the free plate is stressed only by the temperature gradient") {
        CHECK(r.structural);
        CHECK(r.max_von_mises > 0.0);
        CHECK(r.max_reaction < 1e-3);
        CHECK(r.min_stress_zz < 0.0);  // hot side wants to expand, held back in z
        CHECK(r.max_stress_zz > 0.0);
      }
      THEN("results are written for ParaView") {
        REQUIRE(r.files.size() == 2);
        for (const auto& f : r.files) CHECK(fs::exists(f));
      }
    }
  }
}
