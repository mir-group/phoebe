cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)

project(testingPhoebe)

# set global variables required to recognize the package GTest
set (GTEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build")
set (GTEST_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build/lib/libgtest.a")
set (GTEST_MAIN_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/build/lib/libgtest_main.a")
set (GTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest-release-1.10.0/googletest/include")

# include_directories(../src)
# include_directories(../include)

# Locate GTest
find_package(GTest REQUIRED)

# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests example.cpp test.cpp)
target_link_libraries(runTests GTest::GTest GTest::Main)
