#include <bladefem/mesh.hpp>

#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace bladefem {

namespace {

// Gmsh element type ids (MSH 2.2).
constexpr int kLine2 = 1;
constexpr int kTriangle3 = 2;
constexpr int kQuad4 = 3;
constexpr int kTriangle6 = 9;
constexpr int kLine3 = 8;

[[noreturn]] void fail(const std::string& what) {
  throw std::runtime_error("read_gmsh: " + what);
}

void expect_end(std::istream& in, const std::string& section) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("$End" + section, 0) == 0) return;
  }
  fail("missing $End" + section);
}

}  // namespace

Mesh read_gmsh(const std::string& path) {
  std::ifstream in(path);
  if (!in) fail("cannot open '" + path + "'");
  try {
    return read_gmsh(in);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(std::string(e.what()) + " [" + path + "]");
  }
}

Mesh read_gmsh(std::istream& in) {
  Mesh m;
  std::unordered_map<long, int> node_index;  // Gmsh ids are 1-based and may have gaps
  bool have_format = false, have_nodes = false, have_elements = false;

  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();

    if (line == "$MeshFormat") {
      std::string version;
      int file_type = 0, data_size = 0;
      in >> version >> file_type >> data_size;
      if (version.rfind("2.", 0) != 0) {
        fail("unsupported MSH version " + version + "; export with `gmsh -format msh22`");
      }
      if (file_type != 0) fail("binary MSH files are not supported; export as ASCII");
      have_format = true;
      expect_end(in, "MeshFormat");
    } else if (line == "$PhysicalNames") {
      int count = 0;
      in >> count;
      for (int k = 0; k < count; ++k) {
        int dim = 0, tag = 0;
        std::string name;
        in >> dim >> tag;
        std::getline(in, name);
        const auto first = name.find('"'), last = name.rfind('"');
        if (first == std::string::npos || last == first) fail("malformed physical name");
        m.names[name.substr(first + 1, last - first - 1)] = tag;
      }
      expect_end(in, "PhysicalNames");
    } else if (line == "$Nodes") {
      long count = 0;
      in >> count;
      m.nodes.reserve(static_cast<std::size_t>(count));
      for (long k = 0; k < count; ++k) {
        long id = 0;
        double x = 0, y = 0, z = 0;
        if (!(in >> id >> x >> y >> z)) fail("truncated $Nodes section");
        node_index[id] = static_cast<int>(m.nodes.size());
        m.nodes.emplace_back(x, y);
      }
      have_nodes = true;
      expect_end(in, "Nodes");
    } else if (line == "$Elements") {
      if (!have_nodes) fail("$Elements before $Nodes");
      long count = 0;
      in >> count;
      for (long k = 0; k < count; ++k) {
        long id = 0;
        int type = 0, ntags = 0;
        if (!(in >> id >> type >> ntags)) fail("truncated $Elements section");
        int physical = 0;
        for (int t = 0; t < ntags; ++t) {
          int value = 0;
          in >> value;
          if (t == 0) physical = value;
        }
        std::string rest;
        std::getline(in, rest);
        std::istringstream nodes_in(rest);
        const auto node = [&]() {
          long gid = 0;
          if (!(nodes_in >> gid)) fail("element " + std::to_string(id) + " has too few nodes");
          const auto it = node_index.find(gid);
          if (it == node_index.end()) fail("element " + std::to_string(id) + " references unknown node");
          return it->second;
        };

        if (type == kTriangle3) {
          const int a = node(), b = node(), c = node();
          m.triangles.push_back({a, b, c});
          m.triangle_tags.push_back(physical);
        } else if (type == kLine2) {
          const int a = node(), b = node();
          m.edges.push_back({{a, b}, physical});
        } else if (type == kQuad4 || type == kTriangle6 || type == kLine3) {
          fail("only linear triangles are supported (element type " + std::to_string(type) +
               "); mesh with `Mesh.ElementOrder = 1` and no recombination");
        }
        // Points and other auxiliary entities are skipped.
      }
      have_elements = true;
      expect_end(in, "Elements");
    }
  }

  if (!have_format) fail("missing $MeshFormat");
  if (!have_elements || m.triangles.empty()) fail("no triangles found");
  m.orient_counter_clockwise();
  return m;
}

}  // namespace bladefem
