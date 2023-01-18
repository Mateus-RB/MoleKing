#Cmake Requirements
cmake_minimum_required(VERSION 3.15)
#project Title and Version
project(MoleKing VERSION 1.4.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )
configure_file(src/MoleKingConfig.hpp.in MoleKing.hpp)
add_definitions(-DBuild_Python=ON)
#find_package(pybind11 REQUIRED)
if(Build_Python)
message(STATUS "Python is built")
#Set project directories
add_subdirectory(src)
add_subdirectory(src/myMath)
add_subdirectory(src/outputProcess)
add_subdirectory(src/berny)
add_subdirectory(src/chemicalUnits)
add_subdirectory(pybind11)

#add_executable(MoleKing src/main.cpp)
pybind11_add_module(MoleKing src/main.cpp)
#Update serash folders
target_include_directories(MoleKing PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_BINARY_DIR}/src/myMath"
                           "${PROJECT_BINARY_DIR}/src/chemicalUnits"
                           "${PROJECT_BINARY_DIR}/src/outputProcess"
                           "${PROJECT_BINARY_DIR}/src/berny"
                           )
SET_TARGET_PROPERTIES (myMath PROPERTIES COMPILE_FLAGS "-fPIC")
SET_TARGET_PROPERTIES (chemicalUnits PROPERTIES COMPILE_FLAGS "-fPIC")
SET_TARGET_PROPERTIES (outputProcess PROPERTIES COMPILE_FLAGS "-fPIC")
SET_TARGET_PROPERTIES (berny PROPERTIES COMPILE_FLAGS "-fPIC")

target_link_libraries(MoleKing 
                        PUBLIC myMath
                        PUBLIC chemicalUnits
                        PUBLIC outputProcess
                        PUBLIC berny
                        )

target_compile_definitions(MoleKing
                                PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
else()
message(STATUS "Python is not built")
#Set project directories
add_subdirectory(src)
add_subdirectory(src/myMath)
add_subdirectory(src/outputProcess)
add_subdirectory(src/berny)
add_subdirectory(src/chemicalUnits)

add_executable(MoleKing src/main2.cpp)
#pybind11_add_module(MoleKing src/main.cpp)
#Update serash folders
target_include_directories(MoleKing PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_BINARY_DIR}/src/myMath"
                           "${PROJECT_BINARY_DIR}/src/chemicalUnits"
                           "${PROJECT_BINARY_DIR}/src/outputProcess"
                           "${PROJECT_BINARY_DIR}/src/berny"
                           )

target_link_libraries(MoleKing 
                        PUBLIC myMath
                        PUBLIC chemicalUnits
                        PUBLIC outputProcess
                        PUBLIC berny
                        )

endif()