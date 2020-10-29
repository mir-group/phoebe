add_subdirectory(lib/kokkos)

include(ExternalProject)
include(FetchContent)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest
    UPDATE_COMMAND ""
)
FetchContent_MakeAvailable(googletest)

# grabs a copy of nlohmann/json from a lighter repository
# which tracks the original, minus the large test data set. 
# This is recommended by the readme of nlohmann/json
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent.git
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/json_src
    UPDATE_COMMAND ""
)
FetchContent_MakeAvailable(json)

ExternalProject_Add(spglib_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/spglib_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/spglib_src
    GIT_REPOSITORY https://github.com/spglib/spglib.git
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)
ExternalProject_Add(pugixml_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/pugixml_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/pugixml_src
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)
ExternalProject_Add(eigen_dep
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen_build
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen_src
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)

link_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/pugixml_build
    ${CMAKE_CURRENT_BINARY_DIR}/spglib_build
    ${CMAKE_CURRENT_BINARY_DIR}/scalapack_build/lib
)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/pugixml_src/src
    ${CMAKE_CURRENT_BINARY_DIR}/spglib_src/src
    ${CMAKE_CURRENT_BINARY_DIR}/eigen_src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest/googletest/include
    ${CMAKE_CURRENT_BINARY_DIR}/json_src/include
)

