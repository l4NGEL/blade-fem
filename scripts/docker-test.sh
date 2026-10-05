#!/usr/bin/env bash
# Configure, build and test inside the blade-fem:dev image.
# Usage: scripts/docker-test.sh [ctest args...]
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
docker run --rm -v "${root}:/src:ro" -v bladefem-build:/build blade-fem:dev bash -c "
  cmake -S /src -B /build -G Ninja -DBLADEFEM_WERROR=ON > /dev/null &&
  cmake --build /build &&
  cd /build && ctest --output-on-failure $*"
