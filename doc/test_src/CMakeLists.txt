cmake_minimum_required( VERSION 2.8 )

project(PennaGenomeAnimalTesting)

set(CMAKE_CXX_FLAGS "-std=c++0x")

add_library(Penna STATIC genome.cpp animal.cpp)

add_executable(genome-test genome-test.cpp)
add_executable(animal-test animal-test.cpp)

target_link_libraries(genome-test Penna)
target_link_libraries(animal-test Penna)

install(TARGETS genome-test animal-test DESTINATION bin)

install(TARGETS Penna
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES genome.hpp animal.hpp DESTINATION include)


