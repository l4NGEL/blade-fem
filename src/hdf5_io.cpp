#include <bladefem/hdf5_io.hpp>

#include <hdf5.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace bladefem {

namespace fs = std::filesystem;

namespace {

// Owns one HDF5 identifier and releases it with the matching close function.
class Handle {
 public:
  Handle(hid_t id, herr_t (*close)(hid_t), const std::string& what) : id_(id), close_(close) {
    if (id_ < 0) throw std::runtime_error("HDF5: " + what);
  }
  ~Handle() { close_(id_); }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  operator hid_t() const { return id_; }

 private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

void silence_hdf5_errors() {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
}

void write_dataset(hid_t file, const std::string& path, hid_t file_type, hid_t mem_type, const void* data,
                   std::size_t rows, std::size_t cols) {
  const hsize_t dims[2] = {rows, cols};
  Handle space(H5Screate_simple(cols == 1 ? 1 : 2, dims, nullptr), H5Sclose, "dataspace for " + path);
  Handle lcpl(H5Pcreate(H5P_LINK_CREATE), H5Pclose, "link property list");
  H5Pset_create_intermediate_group(lcpl, 1);
  Handle dset(H5Dcreate2(file, path.c_str(), file_type, space, lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Dclose,
              "create " + path);
  if (H5Dwrite(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
    throw std::runtime_error("HDF5: write " + path);
  }
}

std::string dims_string(std::size_t rows, std::size_t cols) {
  return cols == 1 ? std::to_string(rows) : std::to_string(rows) + " " + std::to_string(cols);
}

}  // namespace

void write_results(const std::string& basename, const Mesh& mesh, const std::vector<Field>& fields) {
  const std::size_t N = mesh.node_count(), M = mesh.triangle_count();
  for (const auto& f : fields) {
    if (f.components != 1 && f.components != 2) {
      throw std::invalid_argument("write_results: field '" + f.name + "' must have 1 or 2 components");
    }
    const std::size_t rows = f.location == Field::Location::Node ? N : M;
    if (f.values.size() != rows * static_cast<std::size_t>(f.components)) {
      throw std::invalid_argument("write_results: field '" + f.name + "' has " + std::to_string(f.values.size()) +
                                  " values, expected " + std::to_string(rows * static_cast<std::size_t>(f.components)));
    }
  }

  silence_hdf5_errors();
  const fs::path base(basename);
  if (base.has_parent_path()) fs::create_directories(base.parent_path());
  const std::string h5_path = basename + ".h5";
  const std::string h5_name = base.filename().string() + ".h5";

  {
    Handle file(H5Fcreate(h5_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "create " + h5_path);

    std::vector<double> xy(2 * N);
    for (std::size_t i = 0; i < N; ++i) {
      xy[2 * i] = mesh.nodes[i].x();
      xy[2 * i + 1] = mesh.nodes[i].y();
    }
    write_dataset(file, "/mesh/nodes", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, xy.data(), N, 2);

    std::vector<int> tri(3 * M);
    for (std::size_t t = 0; t < M; ++t) {
      for (std::size_t k = 0; k < 3; ++k) tri[3 * t + k] = mesh.triangles[t][k];
    }
    write_dataset(file, "/mesh/triangles", H5T_STD_I32LE, H5T_NATIVE_INT, tri.data(), M, 3);

    for (const auto& f : fields) {
      const std::size_t rows = f.location == Field::Location::Node ? N : M;
      if (f.components == 1) {
        write_dataset(file, "/fields/" + f.name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, f.values.data(), rows, 1);
      } else {
        // XDMF/ParaView vectors are three-dimensional: pad with z = 0.
        std::vector<double> v3(3 * rows, 0.0);
        for (std::size_t r = 0; r < rows; ++r) {
          v3[3 * r] = f.values[2 * r];
          v3[3 * r + 1] = f.values[2 * r + 1];
        }
        write_dataset(file, "/fields/" + f.name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, v3.data(), rows, 3);
      }
    }
  }

  std::ofstream x(basename + ".xdmf");
  if (!x) throw std::runtime_error("write_results: cannot write " + basename + ".xdmf");
  x << "<?xml version=\"1.0\" ?>\n"
    << "<Xdmf Version=\"3.0\">\n  <Domain>\n    <Grid Name=\"mesh\" GridType=\"Uniform\">\n"
    << "      <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << M << "\">\n"
    << "        <DataItem Dimensions=\"" << M << " 3\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">" << h5_name
    << ":/mesh/triangles</DataItem>\n      </Topology>\n"
    << "      <Geometry GeometryType=\"XY\">\n"
    << "        <DataItem Dimensions=\"" << N << " 2\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << h5_name
    << ":/mesh/nodes</DataItem>\n      </Geometry>\n";
  for (const auto& f : fields) {
    const bool node = f.location == Field::Location::Node;
    const std::size_t rows = node ? N : M;
    x << "      <Attribute Name=\"" << f.name << "\" AttributeType=\"" << (f.components == 1 ? "Scalar" : "Vector")
      << "\" Center=\"" << (node ? "Node" : "Cell") << "\">\n"
      << "        <DataItem Dimensions=\"" << dims_string(rows, f.components == 1 ? 1 : 3)
      << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << h5_name << ":/fields/" << f.name
      << "</DataItem>\n      </Attribute>\n";
  }
  x << "    </Grid>\n  </Domain>\n</Xdmf>\n";
}

Dataset read_dataset(const std::string& h5_path, const std::string& dataset_path) {
  silence_hdf5_errors();
  if (!fs::exists(h5_path)) throw std::runtime_error("read_dataset: no such file " + h5_path);
  Handle file(H5Fopen(h5_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "open " + h5_path);
  Handle dset(H5Dopen2(file, dataset_path.c_str(), H5P_DEFAULT), H5Dclose, "open dataset " + dataset_path);
  Handle space(H5Dget_space(dset), H5Sclose, "dataspace of " + dataset_path);

  const int rank = H5Sget_simple_extent_ndims(space);
  std::vector<hsize_t> dims(static_cast<std::size_t>(rank));
  H5Sget_simple_extent_dims(space, dims.data(), nullptr);

  Dataset out;
  std::size_t total = 1;
  for (hsize_t d : dims) {
    out.shape.push_back(static_cast<std::size_t>(d));
    total *= static_cast<std::size_t>(d);
  }
  out.values.resize(total);
  if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.values.data()) < 0) {
    throw std::runtime_error("HDF5: read " + dataset_path);
  }
  return out;
}

}  // namespace bladefem
