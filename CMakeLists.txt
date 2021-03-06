if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(Kokkos_ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "don't touch")
    set(CMAKE_CXX_COMPILER ${CMAKE_CURRENT_SOURCE_DIR}/lib/kokkos/bin/nvcc_wrapper)
endif()

cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 14)

project(Phoebe C CXX Fortran)

## load custom cmake find_package
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -ggdb -fno-omit-frame-pointer -Wno-deprecated-copy -Wno-ignored-attributes")
    add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
endif()

############ OPTIONS ################
# these can be toggled in a build script by 
# supplying arguments to cmake like cmake -DHDF5_AVAIL=OFF ../
option(MPI_AVAIL "Build with MPI wrappers" ON)
option(OMP_AVAIL "Build with OMP" ON)
option(HDF5_AVAIL "Build with HDF5" ON)
option(HDF5_SERIAL "Force build to accomodate serial only HDF5, but still use MPI" OFF)
option(BUILD_DOC "Build documentation" ON)

############### SOURCE ###############

FILE(GLOB_RECURSE SOURCE_FILES src src/*.cpp)
list(REMOVE_ITEM SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

message(${CMAKE_CURRENT_SOURCE_DIR})
message(${SOURCE_FILES})

FILE(GLOB TEST_SOURCES test/*.cpp test/interaction3ph/*.cpp)

include_directories(src src/algebra src/apps src/bands src/bte src/constants src/harmonic src/interaction src/mpi src/observable src/parser src/points)

include(lib/CMakeLists.txt)
include(GoogleTest)

add_executable(phoebe src/main.cpp ${SOURCE_FILES})
add_executable(runTests ${TEST_SOURCES} ${SOURCE_FILES})
set_target_properties(runTests PROPERTIES EXCLUDE_FROM_ALL TRUE)

# dependencies
add_dependencies(phoebe spglib_dep pugixml_dep eigen_dep)
add_dependencies(runTests spglib_dep pugixml_dep eigen_dep)

target_link_libraries(phoebe symspg pugixml kokkos nlohmann_json::nlohmann_json)
enable_testing()
target_link_libraries(runTests symspg pugixml gtest_main kokkos nlohmann_json::nlohmann_json)

gtest_discover_tests(
    runTests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

################ PARALLEL ####################
# build with MPI and scalapack
if (MPI_AVAIL)
  find_package(MPI REQUIRED)
  target_link_libraries(phoebe MPI::MPI_CXX)
  target_link_libraries(runTests MPI::MPI_CXX)
  add_definitions("-DMPI_AVAIL") # this sets turns on ifdefs for mpi in the cpp

  find_library(SCALAPACK_LIB NAMES scalapack scalapack-openmpi mkl_scalapack_lp64 PATHS ENV LD_LIBRARY_PATH)
  if(${SCALAPACK_LIB} MATCHES NOTFOUND)
      ExternalProject_Add(scalapack_dep
          BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/scalapack_build
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/scalapack_src
          GIT_REPOSITORY https://github.com/Reference-ScaLAPACK/scalapack.git
          INSTALL_COMMAND ""
          UPDATE_COMMAND ""
          CMAKE_ARGS
              -DCMAKE_BUILD_TYPE=Release
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
          )
      set(SCALAPACK_LIB scalapack)
      message(WARNING "No SCALAPACK library found in CMAKE_LIBRARY_PATH, building our own.")
      add_dependencies(phoebe scalapack_dep)
  else()
      if(${SCALAPACK_LIB} MATCHES mkl)
          if("${MPI_CXX_LIBRARIES}" MATCHES openmpi)
              find_library(BLACS_LIB NAMES mkl_blacs_openmpi_lp64 PATHS ENV LD_LIBRARY_PATH)
          elseif("${MPI_CXX_LIBRARIES}" MATCHES intel)
              find_library(BLACS_LIB NAMES mkl_blacs_intelmpi_lp64 PATHS ENV LD_LIBRARY_PATH)
          else()
              message(FATAL_ERROR "Confused by MPI library when looking for BLACS.")
          endif()
          if(${BLACS_LIB} MATCHES NOTFOUND)
              message(FATAL_ERROR "Found Intel SCALAPACK but not BLACS")
          endif()
          message(STATUS "Found BLACS: ${BLACS_LIB}")
      endif()
      message(STATUS "Found SCALAPACK: ${SCALAPACK_LIB}")
  endif()

  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  include_directories(${BLAS_INCLUDE_DIR})
  include_directories(${LAPACK_INCLUDE_DIR})
  ## Very important to link scalapack before blas and lapack
  target_link_libraries(phoebe ${SCALAPACK_LIB} ${BLACS_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)
  target_link_libraries(runTests ${SCALAPACK_LIB} ${BLACS_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)

else()

  ## even for the serial version, we require lapack
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  include_directories(${BLAS_INCLUDE_DIR})
  include_directories(${LAPACK_INCLUDE_DIR})
  target_link_libraries(phoebe ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)
  target_link_libraries(runTests ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} -lgfortran)

endif()

# build with openMP
if(OMP_AVAIL)
  find_package(OpenMP)
  target_link_libraries(phoebe OpenMP::OpenMP_CXX)
  target_link_libraries(runTests OpenMP::OpenMP_CXX)
  add_definitions("-DOMP_AVAIL")
endif()

# build with HDF5
if(HDF5_AVAIL)
  add_dependencies(phoebe highfive_dep)
  add_dependencies(runTests highfive_dep)
  target_link_libraries(phoebe hdf5)
  target_link_libraries(runTests hdf5)
  add_definitions("-DHDF5_AVAIL")
  if(HDF5_SERIAL)
    add_definitions("-DHDF5_SERIAL")
  endif()
endif()

############# DOCS ############

find_package(Doxygen)
find_package(Sphinx)

IF (BUILD_DOC)
    IF (NOT DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    else()
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile.out)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen configured")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target(
            docDevel
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc/doxygen
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)

	IF (NOT SPHINX_FOUND)
	  message("Documentation needs also sphinx")
	else()		
	  add_custom_target(
	    doc
	    COMMAND make html
	    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../doc/sphinx
	    COMMENT "Building user documentation with Sphinx"
	    VERBATIM)
	  add_dependencies(doc docDevel)
	endif()

    ENDIF()
endif()



