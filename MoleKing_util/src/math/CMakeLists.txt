add_definitions(-D_USE_MATH_DEFINES)

add_library(math 
Geometry.hpp
Geometry.cpp
MassCenter.hpp
MassCenter.cpp
Matrix.hpp
Matrix.cpp
Vectors.hpp
Vectors.cpp) 

target_include_directories(math PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")