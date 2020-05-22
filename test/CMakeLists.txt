cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)

project(testingPhoebe)

# set global variables required to recognize the package GTest
set (GTEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build")
set (GTEST_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build/lib/libgtest.a")
set (GTEST_MAIN_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build/lib/libgtest_main.a")
set (GTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/googletest/include")

include_directories(../src)
include_directories(../include)
include_directories(../include/Eigen)
include_directories(../include/)
include_directories(../lib/pugixml-1.10/src)

# Locate GTest
find_package(GTest REQUIRED)

SET(GCC_COVERAGE_COMPILE_FLAGS "-Wno-deprecated-copy -Wno-ignored-attributes")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

# Link runTests with what we want to test and the GTest and pthread library
# add_executable(runTests example.cpp test.cpp)
add_executable(runTests test.cpp ../src/constants/periodic_table.cpp ../src/algebra/utilities.cpp ../src/exceptions/exceptions.cpp ../src/io.cpp ../src/statistics.cpp ../src/context.cpp ../src/crystal.cpp ../src/harmonic/harmonic.cpp ../src/harmonic/phonon_h0.cpp ../src/points.cpp ../src/state.cpp ../src/active_bandstructure.cpp ../src/harmonic/electron_h0_fourier.cpp ../src/harmonic/electron_h0_wannier.cpp ../src/harmonic/window.cpp ../src/bte/vector_bte.cpp ../src/pugixml.cpp ../src/parser/qe_input_parser.cpp ../src/observable/observable.cpp ../src/bte/scattering.cpp ../src/bte/ph_scattering.cpp ../src/bte/drift.cpp ../src/delta_function/delta_function.cpp ../src/statistics_sweep.cpp ../src/parser/ifc3_parser.cpp ../src/interaction/interaction_3ph.cpp pointsTest.cpp)

target_link_libraries(runTests GTest::GTest GTest::Main)

# link spglib
# LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/spglib-1.14.1)
# LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
TARGET_LINK_LIBRARIES(runTests ${CMAKE_CURRENT_SOURCE_DIR}/../lib/spglib-1.14.1/libsymspg.a)
