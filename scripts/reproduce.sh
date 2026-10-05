#!/usr/bin/env bash
# Rebuild everything behind the README from scratch: build, tests, verification
# studies, blade meshes and runs, mesh-independence study and figures.
# Runs on Ubuntu 24.04 with the packages in the Dockerfile, e.g.
#   docker build -t blade-fem:dev . && docker run --rm -v "$PWD:/src" blade-fem:dev scripts/reproduce.sh
set -euo pipefail
cd "$(dirname "$0")/.."

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure

./build/bladefem-verify --out results/verification

python3 scripts/make_blade_geo.py
python3 scripts/make_blade_geo.py --te-passage
gmsh meshes/blade.geo -2 -format msh22 -o meshes/blade.msh > /dev/null
gmsh meshes/blade_te.geo -2 -format msh22 -o meshes/blade_te.msh > /dev/null
./build/blade-fem cases/blade.json
./build/blade-fem cases/blade_te_passage.json

python3 scripts/blade_mesh_study.py --bin build/blade-fem
python3 scripts/plot_results.py
