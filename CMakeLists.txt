cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)

project(Phoebe)

# load externalProject module
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
# install libraries
ExternalProject_Add(lib
    PREFIX lib
    TMP_DIR lib
    SOURCE_DIR lib
    BINARY_DIR lib
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/lib
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

# set global variables required to recognize the package GTest
set (GTEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0")
set (GTEST_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgtest.a")
set (GTEST_MAIN_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgtest_main.a")
set (GTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/googletest/include")

include_directories(include)
include_directories(lib/pugixml-1.10/build/include)
include_directories(lib/eigen-3.3.7)
include_directories(lib/spglib-1.14.1/build/include)
include_directories(lib/googletest-release-1.10.0)
include_directories(lib/googletest-release-1.10.0/googletest/include)
include_directories(lib/googletest-release-1.10.0/build/include.gtest)

# Locate GTest
find_package(GTest REQUIRED)

SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -O3 -g -Wno-deprecated-copy -Wno-ignored-attributes")
# Note: no-dep-copy and no-ign-attr are to suppress warnings from Eigen library
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

# Link runTests with what we want to test and the GTest and pthread library

SET (SOURCE_FILES src/apps/app.cpp src/apps/phonon_transport_app.cpp src/apps/polarization_app.cpp src/constants/periodic_table.cpp src/algebra/utilities.cpp src/exceptions/exceptions.cpp src/io.cpp src/statistics.cpp src/context.cpp src/crystal.cpp src/harmonic/harmonic.cpp src/harmonic/phonon_h0.cpp src/points.cpp src/state.cpp src/active_bandstructure.cpp src/harmonic/electron_h0_fourier.cpp src/harmonic/electron_h0_wannier.cpp src/harmonic/window.cpp src/bte/vector_bte.cpp src/parser/qe_input_parser.cpp src/observable/observable.cpp src/bte/scattering.cpp src/bte/ph_scattering.cpp src/bte/drift.cpp src/delta_function/delta_function.cpp src/statistics_sweep.cpp src/parser/ifc3_parser.cpp src/interaction/interaction_3ph.cpp src/apps/dos_app.cpp src/apps/bands_app.cpp src/observable/phonon_thermal_cond.cpp src/observable/phonon_viscosity.cpp)
FILE(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/*.cpp ${CMAKE_SOURCE_DIR}/test/interaction3ph/*.cpp)

add_executable(phoebe src/main.cpp ${SOURCE_FILES})
add_executable(runTests ${TEST_SOURCES} ${SOURCE_FILES})

add_dependencies(phoebe lib)
add_dependencies(runTests lib)

# Add Google Test to testing executable
# target_link_libraries(runTests STATIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgtest.a ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgtest_main.a)
# target_link_libraries(runTests GTest::Main)

target_link_libraries(runTests -lpthread ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgmock_main.a ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0/build/lib/libgmock.a)

# link spglib
target_link_libraries(phoebe ${CMAKE_CURRENT_SOURCE_DIR}/lib/spglib-1.14.1/build/lib/libsymspg.a)
target_link_libraries(runTests ${CMAKE_CURRENT_SOURCE_DIR}/lib/spglib-1.14.1/build/lib/libsymspg.a)

# link pugixml
target_link_libraries(phoebe ${CMAKE_CURRENT_SOURCE_DIR}/lib/pugixml-1.10/build/lib/libpugixml.a)
target_link_libraries(runTests ${CMAKE_CURRENT_SOURCE_DIR}/lib/pugixml-1.10/build/lib/libpugixml.a)

############# DOCS ############

# check if Doxygen is installed
find_package(Doxygen)

# we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

IF (BUILD_DOC)
   IF (NOT DOXYGEN_FOUND)
      message("Doxygen need to be installed to generate the doxygen documentation")
   ENDIF()
   # set input and output files
   set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)
   set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile.out)

   # request to configure the file
   configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
   message("Doxygen configured")

   # note the option ALL which allows to build the docs together with the application
   add_custom_target(
	doc
	COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
   	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
