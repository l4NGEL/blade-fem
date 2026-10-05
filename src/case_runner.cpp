#include <bladefem/case_runner.hpp>
#include <bladefem/elastic2d.hpp>
#include <bladefem/hdf5_io.hpp>
#include <bladefem/thermal2d.hpp>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace bladefem {

CaseReport run_case(const CaseConfig& c, const Mesh& mesh) {
  CaseReport r;
  r.nodes = mesh.node_count();
  r.triangles = mesh.triangle_count();
  std::vector<Field> fields;

  Eigen::VectorXd temperature;
  if (c.thermal) {
    ThermalProblem tp;
    tp.conductivity = c.thermal->conductivity;
    if (c.thermal->source != 0.0) {
      const double q = c.thermal->source;
      tp.source = [q](const Vec2&) { return q; };
    }
    for (const auto& f : c.thermal->temperature) {
      const double T = f.T;
      tp.dirichlet.push_back({f.boundary, [T](const Vec2&) { return T; }});
    }
    for (const auto& cv : c.thermal->convection) tp.convection.push_back({cv.boundary, cv.h, cv.T_inf});
    for (const auto& fl : c.thermal->flux) tp.flux.push_back({fl.boundary, fl.q_in});

    const ThermalSolution ts = solve_thermal(mesh, tp);
    temperature = ts.temperature;
    r.thermal = true;
    Eigen::Index imax = 0;
    r.T_max = temperature.maxCoeff(&imax);
    r.T_min = temperature.minCoeff();
    r.T_max_at = mesh.nodes[static_cast<std::size_t>(imax)];
    r.boundary_heat = ts.boundary_heat;
    fields.push_back({"temperature", Field::Location::Node, 1,
                      std::vector<double>(temperature.data(), temperature.data() + temperature.size())});
  }

  if (c.structural) {
    const auto& s = *c.structural;
    ElasticProblem ep;
    ep.material = s.material;
    ep.model = s.model;
    ep.T_ref = s.T_ref;
    ep.temperature = temperature;
    for (const auto& sup : s.supports) ep.supports.push_back({sup.at, sup.x, sup.y});
    for (const auto& d : s.displacement) {
      if (d.ux) {
        const double v = *d.ux;
        ep.displacement.push_back({d.boundary, 0, [v](const Vec2&) { return v; }});
      }
      if (d.uy) {
        const double v = *d.uy;
        ep.displacement.push_back({d.boundary, 1, [v](const Vec2&) { return v; }});
      }
    }
    for (const auto& l : s.traction) {
      const Vec2 t = l.traction;
      ep.traction.push_back({l.boundary, [t](const Vec2&) { return t; }});
    }

    const ElasticSolution es = solve_elastic(mesh, ep);
    r.structural = true;
    const auto vm_max = std::max_element(es.von_mises.begin(), es.von_mises.end());
    r.max_von_mises = *vm_max;
    r.max_von_mises_at = mesh.centroid(static_cast<std::size_t>(vm_max - es.von_mises.begin()));
    r.min_stress_zz = *std::min_element(es.stress_zz.begin(), es.stress_zz.end());
    r.max_stress_zz = *std::max_element(es.stress_zz.begin(), es.stress_zz.end());
    for (std::size_t i = 0; i < mesh.node_count(); ++i) {
      r.max_displacement = std::max(r.max_displacement, es.node_displacement(static_cast<int>(i)).norm());
    }
    r.axial_strain = es.axial_strain;
    r.max_reaction = es.reactions.cwiseAbs().maxCoeff();

    std::vector<double> sxx(r.triangles), syy(r.triangles), sxy(r.triangles);
    for (std::size_t t = 0; t < r.triangles; ++t) {
      sxx[t] = es.stress[t][0];
      syy[t] = es.stress[t][1];
      sxy[t] = es.stress[t][2];
    }
    fields.push_back({"displacement", Field::Location::Node, 2,
                      std::vector<double>(es.displacement.data(), es.displacement.data() + es.displacement.size())});
    fields.push_back({"von_mises", Field::Location::Cell, 1, es.von_mises});
    fields.push_back({"von_mises_nodal", Field::Location::Node, 1, nodal_average(mesh, es.von_mises)});
    fields.push_back({"stress_xx", Field::Location::Cell, 1, sxx});
    fields.push_back({"stress_yy", Field::Location::Cell, 1, syy});
    fields.push_back({"stress_xy", Field::Location::Cell, 1, sxy});
    fields.push_back({"stress_zz", Field::Location::Cell, 1, es.stress_zz});
    fields.push_back({"stress_zz_nodal", Field::Location::Node, 1, nodal_average(mesh, es.stress_zz)});
  }

  if (!c.output.empty()) {
    write_results(c.output, mesh, fields);
    r.files = {c.output + ".h5", c.output + ".xdmf"};
  }
  return r;
}

void print_report(std::ostream& out, const CaseConfig& c, const CaseReport& r) {
  const auto mm = [](const Vec2& p) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << '(' << p.x() * 1e3 << ", " << p.y() * 1e3 << ") mm";
    return s.str();
  };
  const auto saved = out.flags();
  out << (c.name.empty() ? "case" : c.name) << "\n"
      << "  mesh            " << r.nodes << " nodes, " << r.triangles << " triangles\n";
  if (r.thermal) {
    out << std::fixed << std::setprecision(1) << "  temperature     " << r.T_min << " .. " << r.T_max << " K (max at "
        << mm(r.T_max_at) << ")\n";
    double net = 0.0;
    for (const auto& [name, q] : r.boundary_heat) {
      out << "  heat in         " << std::setw(16) << std::left << name << std::right << std::setprecision(2) << q
          << " W/m\n";
      net += q;
    }
    out << std::scientific << std::setprecision(2) << "  heat balance    " << net << " W/m\n";
  }
  if (r.structural) {
    out << std::fixed << std::setprecision(1) << "  model           " << to_string(c.structural->model) << "\n"
        << "  von Mises max   " << r.max_von_mises / 1e6 << " MPa at " << mm(r.max_von_mises_at) << "\n"
        << "  sigma_zz        " << r.min_stress_zz / 1e6 << " .. " << r.max_stress_zz / 1e6 << " MPa\n"
        << std::setprecision(2) << "  |u| max         " << r.max_displacement * 1e6 << " um\n";
    if (c.structural->model == PlaneModel::GeneralizedStrain) {
      out << std::scientific << std::setprecision(3) << "  axial strain    e0 = " << r.axial_strain[0]
          << ", kx = " << r.axial_strain[1] << " 1/m, ky = " << r.axial_strain[2] << " 1/m\n";
    }
    out << std::scientific << std::setprecision(2) << "  max reaction    " << r.max_reaction << " N/m\n";
  }
  for (const auto& f : r.files) out << "  wrote           " << f << "\n";
  out.flags(saved);
}

}  // namespace bladefem
