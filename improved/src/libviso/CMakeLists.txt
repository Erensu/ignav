# project
cmake_minimum_required (VERSION 2.6)
project (libviso)

# directories
set(LIBVISO2_SRC_DIR src/libvisio)

# include directory
include_directories("${LIBVISO2_SRC_DIR}")

# use sse3 instruction set
SET(CMAKE_CXX_FLAGS "-msse3")

# sources
FILE(GLOB LIBVISO2_SRC_FILES "*.cc")

# demo program
add_library(viso ${LIBVISO2_SRC_FILES})
target_link_libraries (viso png)

