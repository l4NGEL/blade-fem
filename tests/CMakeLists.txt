find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
  include(FetchContent)
  FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
include(Catch)

file(GLOB BLADEFEM_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp)
add_executable(bladefem_tests ${BLADEFEM_TEST_SOURCES})
target_link_libraries(bladefem_tests PRIVATE bladefem Catch2::Catch2WithMain)
bladefem_warnings(bladefem_tests)
catch_discover_tests(bladefem_tests)
