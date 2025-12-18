#add_definitions(-D_USE_MATH_DEFINES)
set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )
add_library(myMath 
Geometry.cpp
MassCenter.cpp
Matrix.cpp
Vectors.cpp) 

#target_include_directories(math PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")