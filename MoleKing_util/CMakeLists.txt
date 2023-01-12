#Cmake Requirements
cmake_minimum_required(VERSION 3.15)
#project Title and Version
project(MoleKing_Utils VERSION 1.4.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(src/MoleKing_utilConfig.hpp.in MoleKing_util.hpp)

#Set project directories
#add_subdirectory(pybind11)
add_subdirectory(src)
add_subdirectory(src/myMath)
add_subdirectory(src/outputProcess)
add_subdirectory(src/berny)
add_subdirectory(src/chemicalUnits)

add_executable(MoleKing_util src/main.cpp)

#Update serash folders
target_include_directories(MoleKing_util PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_BINARY_DIR}/src/myMath"
                           "${PROJECT_BINARY_DIR}/src/chemicalUnits"
                           "${PROJECT_BINARY_DIR}/src/outputProcess"
                           "${PROJECT_BINARY_DIR}/src/berny"
                           )

target_link_libraries(MoleKing_util 
                        PUBLIC myMath
                        PUBLIC chemicalUnits
                        PUBLIC outputProcess
                        PUBLIC berny)