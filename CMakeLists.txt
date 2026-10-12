cmake_minimum_required(VERSION 3.20)
project(blade_fem
  VERSION 0.1.0
  DESCRIPTION "Finite element thermal and thermal-stress analysis of turbine-blade sections"
  LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(BLADEFEM_BUILD_TESTS "Build the Catch2 test suite" ${PROJECT_IS_TOP_LEVEL})
option(BLADEFEM_BUILD_APPS "Build the command-line tools" ON)
option(BLADEFEM_WERROR "Treat compiler warnings as errors" OFF)
option(BLADEFEM_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(BLADEFEM_COVERAGE "Instrument for gcov line coverage (GCC/Clang)" OFF)

if(BLADEFEM_SANITIZE AND NOT MSVC)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=address,undefined)
endif()
if(BLADEFEM_COVERAGE AND NOT MSVC)
  add_compile_options(--coverage -O0 -g)
  add_link_options(--coverage)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.71 REQUIRED COMPONENTS program_options)
find_package(HDF5 REQUIRED COMPONENTS C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

function(bladefem_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
  endif()
  if(BLADEFEM_WERROR)
    target_compile_options(${target} PRIVATE $<IF:$<CXX_COMPILER_ID:MSVC>,/WX,-Werror>)
  endif()
endfunction()

# ---- library ---------------------------------------------------------------

configure_file(include/bladefem/version.hpp.in ${PROJECT_BINARY_DIR}/generated/bladefem/version.hpp @ONLY)

file(GLOB BLADEFEM_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(bladefem ${BLADEFEM_SOURCES})
add_library(bladefem::bladefem ALIAS bladefem)
target_compile_features(bladefem PUBLIC cxx_std_17)
target_include_directories(bladefem PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(bladefem PUBLIC Eigen3::Eigen PRIVATE Boost::headers HDF5::HDF5)
set_target_properties(bladefem PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
bladefem_warnings(bladefem)

# ---- command-line tools ------------------------------------------------------

# apps/<name>.cpp -> executable <name> (underscores become dashes).
set(BLADEFEM_APP_TARGETS)
if(BLADEFEM_BUILD_APPS)
  file(GLOB BLADEFEM_APPS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/apps/*.cpp)
  foreach(app_src ${BLADEFEM_APPS})
    get_filename_component(app_name ${app_src} NAME_WE)
    string(REPLACE "_" "-" app_name ${app_name})
    add_executable(${app_name} ${app_src})
    target_link_libraries(${app_name} PRIVATE bladefem Boost::program_options)
    bladefem_warnings(${app_name})
    list(APPEND BLADEFEM_APP_TARGETS ${app_name})
  endforeach()
endif()

# ---- install / package -------------------------------------------------------

install(TARGETS bladefem ${BLADEFEM_APP_TARGETS}
  EXPORT bladefemTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/bladefem DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hpp")
install(FILES ${PROJECT_BINARY_DIR}/generated/bladefem/version.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bladefem)

set(BLADEFEM_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/bladefem)
install(EXPORT bladefemTargets NAMESPACE bladefem:: DESTINATION ${BLADEFEM_CMAKE_DIR})
configure_package_config_file(cmake/bladefemConfig.cmake.in ${PROJECT_BINARY_DIR}/bladefemConfig.cmake
  INSTALL_DESTINATION ${BLADEFEM_CMAKE_DIR})
write_basic_package_version_file(${PROJECT_BINARY_DIR}/bladefemConfigVersion.cmake COMPATIBILITY SameMinorVersion)
install(FILES ${PROJECT_BINARY_DIR}/bladefemConfig.cmake ${PROJECT_BINARY_DIR}/bladefemConfigVersion.cmake
  DESTINATION ${BLADEFEM_CMAKE_DIR})

# ---- tests -------------------------------------------------------------------

if(BLADEFEM_BUILD_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()
