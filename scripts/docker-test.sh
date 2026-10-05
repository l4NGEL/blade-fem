#!/usr/bin/env bash
# Configure, build and test inside the blade-fem:dev image (docker build -t blade-fem:dev .).
# Usage: scripts/docker-test.sh [ctest args...]
set -euo pipefail
cd "$(dirname "$0")/.."
root="$(pwd -W 2>/dev/null || pwd)"   # Windows path under Git Bash, POSIX path elsewhere
export MSYS_NO_PATHCONV=1
docker run --rm -v "${root}:/src:ro" -v bladefem-build:/build blade-fem:dev bash -c "
  cmake -S /src -B /build -G Ninja -DBLADEFEM_WERROR=ON > /dev/null &&
  cmake --build /build &&
  cd /build && ctest --output-on-failure $*"
