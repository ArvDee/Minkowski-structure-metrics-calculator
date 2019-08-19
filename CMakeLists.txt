cmake_minimum_required (VERSION 2.8.12)
# project name 
project (FloppyBox)

# Set debug build type for gdb flags
set(CMAKE_BUILD_TYPE Debug)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CODE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CODE_GJK_DIR ${CMAKE_SOURCE_DIR}/src/gjk)
set(CODE_SHAPE_DIR ${CMAKE_SOURCE_DIR}/src/gjk/shape)
set(CODE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Include all external headers from the "include" directory
include_directories(${CODE_INCLUDE_DIR} ${CODE_GJK_DIR} ${CODE_SHAPE_DIR})

# All sources files are from the "src" directory
file(GLOB SOURCES ${CODE_SOURCE_DIR}/*.cpp )
file(GLOB GJK_SOURCES ${CODE_GJK_DIR}/*.cpp )
file(GLOB SHAPE_SOURCES ${CODE_SHAPE_DIR}/*.cpp )

# Include the Voro++ header and set the library path
include_directories(/usr/local/include/voro++)
link_directories(/usr/local/lib)

# flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Additional arguments
add_compile_options(-Wall -lm -O3 -g)

# add the executable
add_executable(fbmc ${SOURCES} ${GJK_SOURCES} ${SHAPE_SOURCES})

# link math library
target_link_libraries(fbmc m)
# link Voro++ library
target_link_libraries(fbmc voro++)

# Install
install(TARGETS fbmc DESTINATION ${CMAKE_SOURCE_DIR})
