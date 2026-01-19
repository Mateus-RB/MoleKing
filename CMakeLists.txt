#Cmake Requirements
cmake_minimum_required(VERSION 3.15...4.0)
#project Title and Version
project(MoleKing VERSION 1.6.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if (MSVC)
  # Statically link the MSVC runtime to avoid missing DLLs in wheels
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
if (MSVC)
  add_compile_options(/W1)
else()
  add_compile_options(-Wall)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
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
pybind11_add_module(_core src/main.cpp)
# Force the Python extension to land where setuptools expects (passed via PY_EXT_OUTPUT_DIR)
if(PY_EXT_OUTPUT_DIR)
  foreach(_cfg IN ITEMS "" DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    string(TOUPPER "${_cfg}" _ucfg)
    if(_cfg STREQUAL "")
      set_target_properties(_core PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PY_EXT_OUTPUT_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${PY_EXT_OUTPUT_DIR}
        ARCHIVE_OUTPUT_DIRECTORY ${PY_EXT_OUTPUT_DIR}
      )
    else()
      set_target_properties(_core PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_${_ucfg} ${PY_EXT_OUTPUT_DIR}
        LIBRARY_OUTPUT_DIRECTORY_${_ucfg} ${PY_EXT_OUTPUT_DIR}
        ARCHIVE_OUTPUT_DIRECTORY_${_ucfg} ${PY_EXT_OUTPUT_DIR}
      )
    endif()
  endforeach()
endif()
#Update serash folders
target_include_directories(_core PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_BINARY_DIR}/src/myMath"
                           "${PROJECT_BINARY_DIR}/src/chemicalUnits"
                           "${PROJECT_BINARY_DIR}/src/outputProcess"
                           "${PROJECT_BINARY_DIR}/pybind11"
                           )
set_target_properties(myMath chemicalUnits outputProcess berny
  PROPERTIES POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(_core 
                        PUBLIC myMath
                        PUBLIC chemicalUnits
                        PUBLIC outputProcess)

target_compile_definitions(_core
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
                        PUBLIC berny)
                        
enable_testing()
# First Test
add_test(NAME Usage COMMAND MoleKing)
set_tests_properties(Usage 
                    PROPERTIES PASS_REGULAR_EXPRESSION "Usage: .* parameter"
                    )
# Second Test
add_test(NAME Version COMMAND MoleKing --version)
set_tests_properties(Version 
                    PROPERTIES PASS_REGULAR_EXPRESSION "MoleKing.*Version.*${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
                    )
# Third Test
add_test(NAME Water COMMAND MoleKing water)
set_tests_properties(Water 
                    PROPERTIES PASS_REGULAR_EXPRESSION "Molecule.*O.*1.*H.*2.*\n.*0.*\n.*1.*\n.*0.*\n.*2.*"
                    )
# Fourth Test
add_test(NAME Methanol COMMAND MoleKing methanol)
set_tests_properties(Methanol 
                    PROPERTIES PASS_REGULAR_EXPRESSION ".*1.*0.*2.*\n.*1.*0.*3.*\n.*1.*0.*4.*\n.*2.*0.*3.*\n.*2.*0.*4.*\n.*3.*0.*4.*\n.*0.*4.*5.*\n"
                    )
endif()
