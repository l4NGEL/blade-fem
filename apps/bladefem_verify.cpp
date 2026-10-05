// bladefem-verify: mesh-refinement studies against closed-form solutions.
// Writes one CSV per study and prints observed convergence rates.
//
//   bladefem-verify --out results/verification

#include <bladefem/elastic2d.hpp>
#include <bladefem/mesh.hpp>
#include <bladefem/thermal2d.hpp>
#include <bladefem/verification.hpp>

#include <boost/program_options.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace po = boost::program_options;
namespace fs = std::filesystem;
using namespace bladefem;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Row {
  double h;
  std::vector<double> values;
};

// Writes the study to CSV and prints it with the observed order of each error
// column between consecutive levels.
void report(const fs::path& dir, const std::string& name, const std::string& title,
            const std::vector<std::string>& columns, const std::vector<Row>& rows, std::size_t rate_columns) {
  std::ofstream csv(dir / (name + ".csv"));
  csv << "h";
  for (const auto& c : columns) csv << ',' << c;
  csv << '\n' << std::setprecision(10);
  for (const auto& r : rows) {
    csv << r.h;
    for (double v : r.values) csv << ',' << v;
    csv << '\n';
  }

  std::cout << "\n" << title << "\n" << std::setw(12) << "h";
  for (std::size_t c = 0; c < columns.size(); ++c) {
    std::cout << std::setw(16) << columns[c];
    if (c < rate_columns) std::cout << std::setw(7) << "rate";
  }
  std::cout << '\n';
  for (std::size_t i = 0; i < rows.size(); ++i) {
    std::cout << std::scientific << std::setprecision(3) << std::setw(12) << rows[i].h;
    for (std::size_t c = 0; c < columns.size(); ++c) {
      std::cout << std::setw(16) << rows[i].values[c];
      if (c < rate_columns) {
        if (i == 0) {
          std::cout << std::setw(7) << "-";
        } else {
          std::cout << std::fixed << std::setprecision(2) << std::setw(7)
                    << verification::observed_rate(rows[i - 1].values[c], rows[i].values[c], rows[i - 1].h, rows[i].h)
                    << std::scientific << std::setprecision(3);
        }
      }
    }
    std::cout << '\n';
  }
}

std::vector<Row> heat_manufactured() {
  const double k = 4.0;
  const auto exact = [](const Vec2& p) { return std::sin(kPi * p.x()) * std::sin(kPi * p.y()); };
  const auto grad = [](const Vec2& p) {
    return Vec2(kPi * std::cos(kPi * p.x()) * std::sin(kPi * p.y()), kPi * std::sin(kPi * p.x()) * std::cos(kPi * p.y()));
  };
  std::vector<Row> rows;
  for (int n : {8, 16, 32, 64, 128}) {
    const Mesh m = meshgen::rectangle(0.0, 0.0, 1.0, 1.0, n, n);
    ThermalProblem p;
    p.conductivity = k;
    p.source = [&](const Vec2& x) { return 2.0 * kPi * kPi * k * exact(x); };
    for (const char* side : {"bottom", "right", "top", "left"}) p.dirichlet.push_back({side, [](const Vec2&) { return 0.0; }});
    const auto s = solve_thermal(m, p);
    rows.push_back({m.max_edge_length(),
                    {verification::l2_error(m, s.temperature, exact), verification::h1_seminorm_error(m, s.temperature, grad)}});
  }
  return rows;
}

std::vector<Row> cantilever() {
  const verification::TimoshenkoBeam beam{48.0, 12.0, 3.0e7, 0.3, 1000.0};
  const auto ux = [&](const Vec2& p) { return beam.displacement(p).x(); };
  const auto uy = [&](const Vec2& p) { return beam.displacement(p).y(); };
  std::vector<Row> rows;
  for (int nx : {16, 32, 64, 128, 256}) {
    const Mesh m = meshgen::rectangle(0.0, -beam.D / 2, beam.L, beam.D / 2, nx, nx / 4);
    ElasticProblem p;
    p.material = {beam.E, beam.nu, 0.0};
    p.displacement.push_back({"left", 0, ux});
    p.displacement.push_back({"left", 1, uy});
    p.traction.push_back({"right", [&](const Vec2& x) {
                            const auto s = beam.stress(x);
                            return Vec2(s[0], s[2]);
                          }});
    const auto s = solve_elastic(m, p);
    Eigen::VectorXd vx(m.node_count()), vy(m.node_count());
    for (std::size_t i = 0; i < m.node_count(); ++i) {
      vx[static_cast<Eigen::Index>(i)] = s.displacement[static_cast<Eigen::Index>(2 * i)];
      vy[static_cast<Eigen::Index>(i)] = s.displacement[static_cast<Eigen::Index>(2 * i + 1)];
    }
    const double err = std::hypot(verification::l2_error(m, vx, ux), verification::l2_error(m, vy, uy)) /
                       std::hypot(verification::l2_norm(m, ux), verification::l2_norm(m, uy));
    const double tip = s.node_displacement(m.nearest_node(Vec2(beam.L, 0.0))).y();
    rows.push_back({m.max_edge_length(), {err, std::abs(tip / beam.tip_deflection() - 1.0)}});
  }
  return rows;
}

std::vector<Row> kirsch() {
  const verification::KirschPlate plate{1.0, 1.0};
  std::vector<Row> rows;
  for (int n : {4, 8, 16, 32, 64}) {
    const Mesh m = meshgen::quarter_plate_with_hole(1.0, 5.0, n, 3 * n / 2, 2.0);
    ElasticProblem p;
    p.material = {1000.0, 0.3, 0.0};
    p.displacement.push_back({"symmetry_x", 0, [](const Vec2&) { return 0.0; }});
    p.displacement.push_back({"symmetry_y", 1, [](const Vec2&) { return 0.0; }});
    p.traction.push_back({"right", [&](const Vec2& x) {
                            const auto s = plate.stress(x);
                            return Vec2(s[0], s[2]);
                          }});
    p.traction.push_back({"top", [&](const Vec2& x) {
                            const auto s = plate.stress(x);
                            return Vec2(s[2], s[1]);
                          }});
    const auto s = solve_elastic(m, p);
    double num = 0.0, den = 0.0;
    std::vector<double> sxx(m.triangle_count());
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
      const auto ex = plate.stress(m.centroid(t));
      num += (s.stress[t] - ex).squaredNorm() * m.area(t);
      den += ex.squaredNorm() * m.area(t);
      sxx[t] = s.stress[t][0];
    }
    const double scf = nodal_average(m, sxx)[static_cast<std::size_t>(m.nearest_node(Vec2(0.0, 1.0)))];
    rows.push_back({m.max_edge_length(), {std::sqrt(num / den), std::abs(scf / 3.0 - 1.0), scf}});
  }
  return rows;
}

std::vector<Row> cylinder_gps() {
  const double a = 0.010, b = 0.025, k = 20.0, h = 800.0, Ta = 1100.0, Tinf = 400.0;
  const double E = 200e9, nu = 0.3, alpha = 13e-6;
  const double C = h * (Tinf - Ta) / (k / b + h * std::log(b / a));
  const verification::ThermoelasticRing ring{a, b, E, alpha, Ta, C};
  std::vector<Row> rows;
  for (int f : {1, 2, 4, 8, 16}) {
    const Mesh m = meshgen::annulus(a, b, 48 * f, 6 * f);
    ThermalProblem tp;
    tp.conductivity = k;
    tp.dirichlet.push_back({"inner", [&](const Vec2&) { return Ta; }});
    tp.convection.push_back({"outer", h, Tinf});
    const auto ts = solve_thermal(m, tp);

    ElasticProblem ep;
    ep.material = {E, nu, alpha};
    ep.model = PlaneModel::GeneralizedStrain;
    ep.temperature = ts.temperature;
    ep.T_ref = 293.15;
    ep.supports.push_back({Vec2(a, 0.0), true, true});
    ep.supports.push_back({Vec2(-a, 0.0), false, true});
    const auto es = solve_elastic(m, ep);

    double num = 0.0, den = 0.0;
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
      const Eigen::Vector4d ex = ring.stress_long_cylinder(m.centroid(t), nu);
      const Eigen::Vector4d fe(es.stress[t][0], es.stress[t][1], es.stress[t][2], es.stress_zz[t]);
      num += (fe - ex).squaredNorm() * m.area(t);
      den += ex.squaredNorm() * m.area(t);
    }
    const double szz_inner = nodal_average(m, es.stress_zz)[static_cast<std::size_t>(m.nearest_node(Vec2(a, 0.0)))];
    const double exact_inner = ring.stress_long_cylinder(Vec2(a, 0.0), nu)[3];
    rows.push_back({m.max_edge_length(), {std::sqrt(num / den), std::abs(szz_inner / exact_inner - 1.0), szz_inner / 1e6}});
  }
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  po::options_description opts("Options");
  opts.add_options()("help,h", "show this help")(
      "out,o", po::value<std::string>()->default_value("results/verification"), "directory for the CSV files");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "bladefem-verify: " << e.what() << "\n";
    return 2;
  }
  if (vm.count("help")) {
    std::cout << "Usage: bladefem-verify [--out DIR]\n\n" << opts << "\n";
    return 0;
  }
  const fs::path dir = vm["out"].as<std::string>();
  fs::create_directories(dir);

  report(dir, "heat_manufactured", "Heat conduction, manufactured solution sin(pi x) sin(pi y)",
         {"L2_error", "H1_error"}, heat_manufactured(), 2);
  report(dir, "cantilever", "Timoshenko cantilever (plane stress)", {"disp_L2_rel", "tip_rel_err"}, cantilever(), 2);
  report(dir, "kirsch", "Kirsch plate with hole (plane stress)", {"stress_L2_rel", "scf_rel_err", "scf"}, kirsch(), 2);
  report(dir, "cylinder_gps", "Heated long cylinder, coupled conduction + generalized plane strain",
         {"stress_L2_rel", "szz_inner_err", "szz_inner_MPa"}, cylinder_gps(), 2);
  std::cout << "\nCSV files written to " << dir.string() << "\n";
  return 0;
}
