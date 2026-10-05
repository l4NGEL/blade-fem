#pragma once

// Analysis case read from a JSON file (Boost.PropertyTree). Relative paths in the
// file are resolved against the directory that contains it. Example: cases/blade.json.

#include <bladefem/elastic2d.hpp>

#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

namespace bladefem {

struct CaseConfig {
  std::string name;
  std::string mesh;    // Gmsh 2.2 mesh path
  std::string output;  // result basename: <output>.h5 and <output>.xdmf; empty = no files

  struct Thermal {
    double conductivity = 0.0;  // [W/(m K)]
    double source = 0.0;        // uniform volumetric heat [W/m^3]
    struct Fixed { std::string boundary; double T; };
    struct Convection { std::string boundary; double h, T_inf; };
    struct Flux { std::string boundary; double q_in; };
    std::vector<Fixed> temperature;
    std::vector<Convection> convection;
    std::vector<Flux> flux;
  };

  struct Structural {
    ElasticMaterial material;
    PlaneModel model = PlaneModel::Stress;
    double T_ref = 293.15;  // strain-free temperature [K]
    struct Support { Vec2 at; bool x = true, y = true; };
    struct Fixed { std::string boundary; std::optional<double> ux, uy; };
    struct Load { std::string boundary; Vec2 traction; };
    std::vector<Support> supports;
    std::vector<Fixed> displacement;
    std::vector<Load> traction;
  };

  std::optional<Thermal> thermal;
  std::optional<Structural> structural;
};

// Both throw std::runtime_error naming the offending key.
CaseConfig load_case(const std::string& json_path);
CaseConfig parse_case(std::istream& json, const std::string& base_dir);

// "plane_stress" | "plane_strain" | "generalized_plane_strain"
PlaneModel parse_plane_model(const std::string& name);
std::string to_string(PlaneModel model);

}  // namespace bladefem
