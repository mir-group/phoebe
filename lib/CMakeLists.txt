add_subdirectory(lib/kokkos)

include(ExternalProject)
include(FetchContent)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest
)
FetchContent_MakeAvailable(googletest)

ExternalProject_Add(spglib_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/spglib_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/spglib_src
    GIT_REPOSITORY https://github.com/spglib/spglib.git
    INSTALL_COMMAND ""
)
ExternalProject_Add(pugixml_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/pugixml_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/pugixml_src
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    INSTALL_COMMAND ""
)
ExternalProject_Add(eigen_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen_src
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    INSTALL_COMMAND ""
)

link_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/pugixml_build
    ${CMAKE_CURRENT_BINARY_DIR}/spglib_build
)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/pugixml_src/src
    ${CMAKE_CURRENT_BINARY_DIR}/spglib_src/src
    ${CMAKE_CURRENT_BINARY_DIR}/eigen_src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest/googletest/include
)

