find_package(OpenCV REQUIRED)
add_library(seven-eight-point seven-eight-point.cpp precomp.cpp modelest.cpp)

target_link_libraries(seven-eight-point ${OpenCV_LIBS})
