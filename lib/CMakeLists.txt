cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)

project(PhoebeLibs)

################################
# First, download the packages #
################################

set(EIGEN_ZIPFN ${CMAKE_CURRENT_BINARY_DIR}/eigen.tar.gz)
set(EIGEN_OUTFN ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.3.7)
if(NOT EXISTS ${EIGEN_OUTFN})
  if(NOT EXISTS ${EIGEN_ZIPFN})
    file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz ${EIGEN_ZIPFN}
         SHOW_PROGRESS
         EXPECTED_HASH MD5=9e30f67e8531477de4117506fe44669b)
  endif(NOT EXISTS ${EIGEN_ZIPFN})
  # this is OS-agnostic
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${EIGEN_ZIPFN}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif(NOT EXISTS ${EIGEN_OUTFN})

set(TEST_ZIPFN ${CMAKE_CURRENT_BINARY_DIR}/googletest.tar.gz)
set(TEST_OUTFN ${CMAKE_CURRENT_BINARY_DIR}/googletest-release-1.10.0)
if(NOT EXISTS ${TEST_OUTFN})
  if(NOT EXISTS ${TEST_ZIPFN})
    file(DOWNLOAD https://github.com/google/googletest/archive/release-1.10.0.tar.gz ${TEST_ZIPFN}
         SHOW_PROGRESS
         EXPECTED_HASH MD5=ecd1fa65e7de707cd5c00bdac56022cd)
  endif(NOT EXISTS ${TEST_ZIPFN})
  # this is OS-agnostic
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${TEST_ZIPFN}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif(NOT EXISTS ${TEST_OUTFN})

set(PUGI_ZIPFN ${CMAKE_CURRENT_BINARY_DIR}/pugixml.tar.gz)
set(PUGI_OUTFN ${CMAKE_CURRENT_BINARY_DIR}/pugixml-1.10)
if(NOT EXISTS ${PUGI_OUTFN})
  if(NOT EXISTS ${PUGI_ZIPFN})
    file(DOWNLOAD http://github.com/zeux/pugixml/releases/download/v1.10/pugixml-1.10.tar.gz ${PUGI_ZIPFN}
         SHOW_PROGRESS
         EXPECTED_HASH MD5=f97237e9908201c6d8536210747b66af)
  endif(NOT EXISTS ${PUGI_ZIPFN})
  # this is OS-agnostic
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${PUGI_ZIPFN}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif(NOT EXISTS ${PUGI_OUTFN})

set(SPG_ZIPFN ${CMAKE_CURRENT_BINARY_DIR}/spglib.tar.gz)
set(SPG_OUTFN ${CMAKE_CURRENT_BINARY_DIR}/spglib-1.14.1)
if(NOT EXISTS ${SPG_OUTFN})
  if(NOT EXISTS ${SPG_ZIPFN})
    file(DOWNLOAD https://github.com/spglib/spglib/archive/v1.14.1.zip ${SPG_ZIPFN}
         SHOW_PROGRESS
         EXPECTED_HASH MD5=74abe802f664d8056c0e5aebfe80905d)
  endif(NOT EXISTS ${SPG_ZIPFN})
  # this is OS-agnostic
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${SPG_ZIPFN}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif(NOT EXISTS ${SPG_OUTFN})

###################################################
# Eigen and GoogleTest are header-only libraries. #
# Instead, we must compile pugixml and spglib     #
###################################################

# load externalProject module
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
# install Spglib
ExternalProject_Add(spglib
    PREFIX spglib-1.14.1
    TMP_DIR spglib-1.14.1
    SOURCE_DIR spglib-1.14.1
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/spglib-1.14.1/build
    BUILD_COMMAND make
)

# install Pugixml
ExternalProject_Add(pugixml
    PREFIX pugixml-1.10
    TMP_DIR pugixml-1.10
    SOURCE_DIR pugixml-1.10
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/pugixml-1.10/build
    BUILD_COMMAND make
)

# install GoogleTest
ExternalProject_Add(googletest
    PREFIX googletest-release-1.10.0
    TMP_DIR googletest-release-1.10.0
    SOURCE_DIR googletest-release-1.10.0
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/googletest-release-1.10.0/build
    BUILD_COMMAND make
)
