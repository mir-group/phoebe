cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)

project(Phoebe CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -ggdb -fno-omit-frame-pointer -Wno-deprecated-copy -Wno-ignored-attributes")
    add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
endif()

############### SOURCE ###############

SET (SOURCE_FILES src/bands/bandstructure.cpp src/apps/app.cpp src/apps/phonon_transport_app.cpp src/apps/polarization_app.cpp src/constants/periodic_table.cpp src/utilities.cpp src/exceptions/exceptions.cpp src/io.cpp src/particle.cpp src/context.cpp src/crystal.cpp src/harmonic/harmonic.cpp src/harmonic/phonon_h0.cpp src/points.cpp src/state.cpp src/bands/active_bandstructure.cpp src/harmonic/electron_h0_fourier.cpp src/harmonic/electron_h0_wannier.cpp src/bands/window.cpp src/bte/vector_bte.cpp src/bte/basevector_bte.cpp src/parser/qe_input_parser.cpp src/observable/observable.cpp src/bte/scattering.cpp src/bte/ph_scattering.cpp src/bte/drift.cpp src/delta_function.cpp src/statistics_sweep.cpp src/parser/ifc3_parser.cpp src/interaction/interaction_3ph.cpp src/apps/dos_app.cpp src/apps/bands_app.cpp src/observable/phonon_thermal_cond.cpp src/observable/phonon_viscosity.cpp src/observable/specific_heat.cpp src/mpi/mpiHelper.cpp src/mpi/mpiController.cpp src/algebra/PMatrix.cpp src/algebra/Matrix.cpp src/observable/wigner_phonon_thermal_cond.cpp src/observable/onsager.cpp src/interaction/interaction_epa.cpp src/parser/epa_parser.cpp src/apps/transport_epa_app.cpp)
FILE(GLOB TEST_SOURCES test/*.cpp test/interaction3ph/*.cpp)

include_directories(include)

include(lib/CMakeLists.txt)
include(GoogleTest)

add_executable(phoebe src/main.cpp ${SOURCE_FILES})
add_executable(runTests ${TEST_SOURCES} ${SOURCE_FILES})

add_dependencies(phoebe spglib_dep pugixml_dep eigen_dep)
add_dependencies(runTests spglib_dep pugixml_dep eigen_dep)

target_link_libraries(phoebe symspg pugixml)

enable_testing()
target_link_libraries(runTests symspg pugixml gtest_main)
gtest_discover_tests(
    runTests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

################ PARALLEL ####################
option(MPI_AVAIL "Build with MPI wrappers" ON)
if (MPI_AVAIL)
  find_package(MPI REQUIRED)
  target_link_libraries(phoebe MPI::MPI_CXX)
  target_link_libraries(runTests MPI::MPI_CXX)
  add_definitions("-DMPI_AVAIL") # this sets turns on ifdefs for mpi in the cpp

  add_dependencies(phoebe scalapack_dep)

  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  include_directories(${BLAS_INCLUDE_DIR})
  include_directories(${LAPACK_INCLUDE_DIR})
  ## Very important to link scalapack before blas and lapack
  target_link_libraries(phoebe scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)
  target_link_libraries(runTests scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)

else()

  ## even for the serial version, we require lapack
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  include_directories(${BLAS_INCLUDE_DIR})
  include_directories(${LAPACK_INCLUDE_DIR})
  target_link_libraries(phoebe ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)
  target_link_libraries(runTests ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)

endif()

option(OMP_AVAIL "Build with OMP" OFF)
if(OMP_AVAIL)
        find_package(OpenMP)
        target_link_libraries(phoebe OpenMP::OpenMP_CXX)
        target_link_libraries(runTests OpenMP::OpenMP_CXX)
        add_definitions("-DOMP_AVAIL")
endif()

############# DOCS ############

# check if Doxygen is installed
find_package(Doxygen)

# we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

IF (BUILD_DOC)
    IF (NOT DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    else()
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
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    ENDIF()
endif()

