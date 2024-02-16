cmake_minimum_required(VERSION 3.16.3)
project(Proyecto)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -DNDEBUG -frecord-gcc-switches -O3")

find_package(OpenMP REQUIRED)

add_executable(paos paos.cpp)
add_executable(psoa psoa.cpp)
target_link_libraries(psoa PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(paos PUBLIC OpenMP::OpenMP_CXX)
