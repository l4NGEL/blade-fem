#include <bladefem/case_config.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace bladefem {

namespace pt = boost::property_tree;
namespace fs = std::filesystem;

namespace {

[[noreturn]] void fail(const std::string& what) { throw std::runtime_error("case: " + what); }

template <class T>
T required(const pt::ptree& node, const std::string& key, const std::string& where) {
  const auto child = node.get_child_optional(key);
  if (!child) fail("missing '" + where + key + "'");
  const auto value = child->get_value_optional<T>();
  if (!value) fail("'" + where + key + "' has an invalid value '" + child->data() + "'");
  return *value;
}

template <class T>
T optional_value(const pt::ptree& node, const std::string& key, const std::string& where, T fallback) {
  return node.get_child_optional(key) ? required<T>(node, key, where) : fallback;
}

Vec2 required_vec2(const pt::ptree& node, const std::string& key, const std::string& where) {
  const auto child = node.get_child_optional(key);
  if (!child) fail("missing '" + where + key + "'");
  std::vector<double> v;
  for (const auto& item : *child) {
    const auto x = item.second.get_value_optional<double>();
    if (!x) fail("'" + where + key + "' must be an array of two numbers");
    v.push_back(*x);
  }
  if (v.size() != 2) fail("'" + where + key + "' must be an array of two numbers");
  return {v[0], v[1]};
}

// Calls fn(item, "section.key[i].") for every element of an optional JSON array.
template <class Fn>
void for_each_item(const pt::ptree& node, const std::string& key, const std::string& where, Fn fn) {
  const auto child = node.get_child_optional(key);
  if (!child) return;
  std::size_t i = 0;
  for (const auto& item : *child) {
    fn(item.second, where + key + "[" + std::to_string(i++) + "].");
  }
}

std::string resolve(const std::string& base_dir, const std::string& path) {
  if (path.empty()) return path;
  const fs::path p(path);
  return (p.is_absolute() ? p : fs::path(base_dir) / p).lexically_normal().string();
}

}  // namespace

PlaneModel parse_plane_model(const std::string& name) {
  if (name == "plane_stress") return PlaneModel::Stress;
  if (name == "plane_strain") return PlaneModel::Strain;
  if (name == "generalized_plane_strain") return PlaneModel::GeneralizedStrain;
  fail("unknown model '" + name + "' (use plane_stress, plane_strain or generalized_plane_strain)");
}

std::string to_string(PlaneModel model) {
  switch (model) {
    case PlaneModel::Stress: return "plane_stress";
    case PlaneModel::Strain: return "plane_strain";
    case PlaneModel::GeneralizedStrain: return "generalized_plane_strain";
  }
  return "?";
}

CaseConfig parse_case(std::istream& json, const std::string& base_dir) {
  pt::ptree root;
  try {
    pt::read_json(json, root);
  } catch (const pt::json_parser_error& e) {
    fail(std::string("invalid JSON: ") + e.what());
  }

  CaseConfig c;
  c.name = root.get<std::string>("name", "");
  c.mesh = resolve(base_dir, required<std::string>(root, "mesh", ""));
  c.output = resolve(base_dir, root.get<std::string>("output", ""));

  if (const auto th = root.get_child_optional("thermal")) {
    const std::string w = "thermal.";
    CaseConfig::Thermal t;
    t.conductivity = required<double>(*th, "conductivity", w);
    t.source = optional_value<double>(*th, "source", w, 0.0);
    for_each_item(*th, "temperature", w, [&](const pt::ptree& n, const std::string& at) {
      t.temperature.push_back({required<std::string>(n, "boundary", at), required<double>(n, "T", at)});
    });
    for_each_item(*th, "convection", w, [&](const pt::ptree& n, const std::string& at) {
      t.convection.push_back({required<std::string>(n, "boundary", at), required<double>(n, "h", at),
                              required<double>(n, "T_inf", at)});
    });
    for_each_item(*th, "flux", w, [&](const pt::ptree& n, const std::string& at) {
      t.flux.push_back({required<std::string>(n, "boundary", at), required<double>(n, "q_in", at)});
    });
    c.thermal = t;
  }

  if (const auto st = root.get_child_optional("structural")) {
    const std::string w = "structural.";
    CaseConfig::Structural s;
    s.material.E = required<double>(*st, "E", w);
    s.material.nu = required<double>(*st, "nu", w);
    s.material.alpha = optional_value<double>(*st, "alpha", w, 0.0);
    s.model = parse_plane_model(st->get<std::string>("model", "plane_stress"));
    s.T_ref = optional_value<double>(*st, "T_ref", w, 293.15);
    for_each_item(*st, "supports", w, [&](const pt::ptree& n, const std::string& at) {
      s.supports.push_back({required_vec2(n, "at", at), optional_value<bool>(n, "x", at, true),
                            optional_value<bool>(n, "y", at, true)});
    });
    for_each_item(*st, "displacement", w, [&](const pt::ptree& n, const std::string& at) {
      CaseConfig::Structural::Fixed f{required<std::string>(n, "boundary", at), std::nullopt, std::nullopt};
      if (n.get_child_optional("ux")) f.ux = required<double>(n, "ux", at);
      if (n.get_child_optional("uy")) f.uy = required<double>(n, "uy", at);
      if (!f.ux && !f.uy) fail("'" + at + "' fixes neither ux nor uy");
      s.displacement.push_back(f);
    });
    for_each_item(*st, "traction", w, [&](const pt::ptree& n, const std::string& at) {
      s.traction.push_back({required<std::string>(n, "boundary", at), required_vec2(n, "t", at)});
    });
    c.structural = s;
  }

  if (!c.thermal && !c.structural) fail("needs a 'thermal' and/or a 'structural' section");
  return c;
}

CaseConfig load_case(const std::string& json_path) {
  std::ifstream in(json_path);
  if (!in) fail("cannot open '" + json_path + "'");
  try {
    return parse_case(in, fs::absolute(json_path).parent_path().string());
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(std::string(e.what()) + " [" + json_path + "]");
  }
}

}  // namespace bladefem
