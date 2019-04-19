find_package(OpenCV REQUIRED)
add_library(one-point one-point.cpp precomp.cpp modelest.cpp)

target_link_libraries(one-point ${OpenCV_LIBS})
