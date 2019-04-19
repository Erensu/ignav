project(carvig)

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -fPIC")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "-msse3")
set(CMAKE_C_FLAGS "-msse3")
set(CMAKE_C_FLAGS "-mavx")

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    SET(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenCV 3.0 QUIET)
if(NOT OpenCV_FOUND)
    find_package(OpenCV 2.4.3 QUIET)
    if(NOT OpenCV_FOUND)
        message(FATAL_ERROR "OpenCV>2.4.3 not found.")
    endif()
endif()
find_package(Eigen3 3.1.0 REQUIRED)
find_package(PythonLibs REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src/ins-gnss/rcv
        ${PROJECT_SOURCE_DIR}/src/libviso
        ${EIGEN3_INCLUDE_DIR}
)
include_directories(${PYTHON_INCLUDE_DIRS})

file(GLOB SRC_FILES "src/ins-gnss/*.cc" "src/ins-gnss/rcv/*.cc")

add_library(carvig ${SRC_FILES})
target_link_libraries(carvig blas lapack z pthread
                      /usr/local/lib/libpng.so
                      jpeg)

set_target_properties(carvig PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
                      ${PROJECT_SOURCE_DIR}/lib)

# add relative pose estimate library
add_subdirectory(src/relative-pose)

# add libvisio
add_subdirectory(src/libviso)

add_executable(carvignav src/ins-gnss/app/igvonav.cc)
target_link_libraries(carvignav carvig pthread
                      ${OpenCV_LIBS}
                      five-point-nister
                      viso)

add_executable(tc-postpos src/ins-gnss/app/tc-postpos.cc)
target_link_libraries(tc-postpos
                      carvig pthread
                      ${OpenCV_LIBS}
                      five-point-nister
                      viso)

