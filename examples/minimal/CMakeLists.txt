# Uses an installed blade-fem as a library:
#   cmake --install build --prefix /opt/bladefem
#   cmake -S examples/minimal -B build-example -DCMAKE_PREFIX_PATH=/opt/bladefem
cmake_minimum_required(VERSION 3.20)
project(bladefem_minimal LANGUAGES C CXX)

find_package(bladefem 0.1 REQUIRED)

add_executable(minimal main.cpp)
target_link_libraries(minimal PRIVATE bladefem::bladefem)
