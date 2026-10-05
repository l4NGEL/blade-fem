// blade-fem: run a thermal / thermo-mechanical case described by a JSON file.
//
//   blade-fem cases/blade.json
//   blade-fem cases/blade.json --mesh meshes/blade_fine.msh --output results/fine/blade

#include <bladefem/case_config.hpp>
#include <bladefem/case_runner.hpp>
#include <bladefem/mesh.hpp>

#include <boost/program_options.hpp>

#include <chrono>
#include <exception>
#include <iostream>

namespace po = boost::program_options;

int main(int argc, char** argv) {
  po::options_description visible("Options");
  visible.add_options()
      ("help,h", "show this help")
      ("mesh,m", po::value<std::string>(), "override the mesh path given in the case file")
      ("output,o", po::value<std::string>(), "override the result basename (<base>.h5 / <base>.xdmf)")
      ("no-output", "solve and report without writing result files");
  po::options_description hidden;
  hidden.add_options()("case", po::value<std::string>(), "case file");
  po::options_description all;
  all.add(visible).add(hidden);
  po::positional_options_description positional;
  positional.add("case", 1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "blade-fem: " << e.what() << "\n";
    return 2;
  }
  if (vm.count("help") || !vm.count("case")) {
    std::cout << "Usage: blade-fem <case.json> [options]\n\n"
                 "Steady conduction followed by a thermal-stress analysis on a Gmsh 2.2 mesh.\n\n"
              << visible << "\n";
    return vm.count("help") ? 0 : 2;
  }

  try {
    bladefem::CaseConfig config = bladefem::load_case(vm["case"].as<std::string>());
    if (vm.count("mesh")) config.mesh = vm["mesh"].as<std::string>();
    if (vm.count("output")) config.output = vm["output"].as<std::string>();
    if (vm.count("no-output")) config.output.clear();

    const auto t0 = std::chrono::steady_clock::now();
    const bladefem::Mesh mesh = bladefem::read_gmsh(config.mesh);
    const bladefem::CaseReport report = bladefem::run_case(config, mesh);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    bladefem::print_report(std::cout, config, report);
    std::cout << "  wall time       " << seconds << " s\n";
  } catch (const std::exception& e) {
    std::cerr << "blade-fem: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
