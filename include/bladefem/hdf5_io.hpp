#pragma once

// Result output: one HDF5 file holding the mesh and all fields, plus an XDMF
// descriptor next to it so the results open directly in ParaView/VisIt.
//
// Layout of <basename>.h5:
//   /mesh/nodes          (N x 2, float64)
//   /mesh/triangles      (M x 3, int32, zero-based)
//   /fields/<name>       (N or M rows; 1 column for scalars, 3 for vectors with z = 0)

#include <bladefem/mesh.hpp>

#include <string>
#include <vector>

namespace bladefem {

struct Field {
  enum class Location { Node, Cell };

  std::string name;
  Location location = Location::Node;
  int components = 1;          // 1 (scalar) or 2 (in-plane vector, written padded to 3)
  std::vector<double> values;  // row-major, rows x components
};

// Writes <basename>.h5 and <basename>.xdmf; parent directories are created.
void write_results(const std::string& basename, const Mesh& mesh, const std::vector<Field>& fields);

struct Dataset {
  std::vector<std::size_t> shape;
  std::vector<double> values;  // converted to double, row-major
};

Dataset read_dataset(const std::string& h5_path, const std::string& dataset_path);

}  // namespace bladefem
