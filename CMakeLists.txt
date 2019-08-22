cmake_minimum_required (VERSION 2.8.12)
project(BOOP_FROM_SNAPSHOTS)

# Set debug build type for gdb flags
set(CMAKE_BUILD_TYPE Debug)

# Set some folders
set(CODE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# All sources files are from the "src" directory
file(GLOB SOURCES ${CODE_SOURCE_DIR}/*.cpp ${CODE_SOURCE_DIR}/*.c)

# Include the Voro++ header and set the library path
include_directories(/usr/local/include/voro++)
link_directories(/usr/local/lib)

# flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Additional arguments
add_compile_options(-Wall -lm -O3 -g)

# Set the executable name
add_executable(boop ${SOURCES})

# link math and Voro++ libraries
target_link_libraries(boop m)
target_link_libraries(boop voro++)

# Install
install(TARGETS boop DESTINATION ${CMAKE_SOURCE_DIR})
