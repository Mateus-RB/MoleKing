# Require modern CMake so we get up-to-date Python and pybind11 helpers.
cmake_minimum_required(VERSION 3.18)

# Declare the project metadata once so all targets share a name and version.
project(MoleKing VERSION 1.6.0 LANGUAGES CXX)
set(PROJECT_DESCRIPTION "Python extension module for theoretical chemistry workflows")

# Expose our custom CMake modules folder (if any helpers live there).
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Cache the frequently referenced directories for readability.
set(MOLEKING_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(MOLEKING_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Allow toggling optional parts of the build from the command line.
if(DEFINED MOLEKING_BUILD_PYTHON_MODULE)
  set(_moleking_build_python_default ${MOLEKING_BUILD_PYTHON_MODULE})
else()
  set(_moleking_build_python_default ON)
endif()
option(Build_Python "Build the MoleKing pybind11 extension module" ${_moleking_build_python_default})

if(DEFINED MOLEKING_BUILD_CLI)
  set(_moleking_build_cli_default ${MOLEKING_BUILD_CLI})
elseif(NOT Build_Python)
  set(_moleking_build_cli_default ON)
else()
  set(_moleking_build_cli_default OFF)
endif()
option(MOLEKING_BUILD_CLI "Build the optional MoleKing CLI application" ${_moleking_build_cli_default})
option(MOLEKING_ENABLE_TESTS "Enable ctest targets (requires the CLI build)" OFF)
option(MOLEKING_ENABLE_IPO "Enable interprocedural optimization / LTO when available" OFF)

# Use a modern standard across compilers and disable vendor extensions.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make all targets position independent, which shared libraries require.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Silence MSVC-specific CRT warnings that are noisy for cross-platform code.
if(MSVC)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# Allow enabling link-time optimization only if the toolchain supports it.
if(MOLEKING_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
  if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
  else()
    message(WARNING "IPO/LTO is not supported by this toolchain: ${ipo_error}")
  endif()
endif()

# Generate the configuration header that exposes version info to C++.
configure_file(
  "${MOLEKING_SRC_DIR}/MoleKingConfig.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/MoleKing.hpp"
  @ONLY
)

# Shared include paths used by multiple libraries/executables.
set(MOLEKING_INTERFACE_INCLUDE_DIRS
  $<BUILD_INTERFACE:${MOLEKING_SRC_DIR}>
  $<BUILD_INTERFACE:${MOLEKING_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)

# Convenience list for private include/search paths.
set(MOLEKING_PRIVATE_INCLUDE_DIRS
  "${MOLEKING_SRC_DIR}"
  "${MOLEKING_INCLUDE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
)

# Small helper to give every target reasonable warning flags.
function(moleking_apply_warnings target_name)
  if(MSVC)
    target_compile_options(${target_name} PRIVATE /W4 /permissive- /EHsc)
  else()
    target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

# Core math utilities reused everywhere else.
add_library(myMath STATIC
  "${MOLEKING_SRC_DIR}/myMath/Geometry.cpp"
  "${MOLEKING_SRC_DIR}/myMath/MassCenter.cpp"
  "${MOLEKING_SRC_DIR}/myMath/Matrix.cpp"
  "${MOLEKING_SRC_DIR}/myMath/Vectors.cpp"
)
target_include_directories(myMath PUBLIC ${MOLEKING_INTERFACE_INCLUDE_DIRS})
moleking_apply_warnings(myMath)

# Atomic/molecular primitives.
add_library(chemicalUnits STATIC
  "${MOLEKING_SRC_DIR}/chemicalUnits/AtomicScale.cpp"
  "${MOLEKING_SRC_DIR}/chemicalUnits/Molecule.cpp"
  "${MOLEKING_SRC_DIR}/chemicalUnits/OPLSff.cpp"
  "${MOLEKING_SRC_DIR}/chemicalUnits/PeriodicTable.cpp"
  "${MOLEKING_SRC_DIR}/chemicalUnits/SupraMolecule.cpp"
  "${MOLEKING_SRC_DIR}/chemicalUnits/pov.cpp"
)
target_include_directories(chemicalUnits PUBLIC ${MOLEKING_INTERFACE_INCLUDE_DIRS})
target_link_libraries(chemicalUnits PUBLIC myMath)
moleking_apply_warnings(chemicalUnits)

# Quantum output parsers and emitters.
add_library(outputProcess STATIC
  "${MOLEKING_SRC_DIR}/outputProcess/G16Process.cpp"
  "${MOLEKING_SRC_DIR}/outputProcess/Psi4Process.cpp"
)
target_include_directories(outputProcess PUBLIC ${MOLEKING_INTERFACE_INCLUDE_DIRS})
target_link_libraries(outputProcess PUBLIC myMath chemicalUnits)
moleking_apply_warnings(outputProcess)

# Geometry optimization helpers.
add_library(berny STATIC "${MOLEKING_SRC_DIR}/berny/Hessian.cpp")
target_include_directories(berny PUBLIC ${MOLEKING_INTERFACE_INCLUDE_DIRS})
target_link_libraries(berny PUBLIC myMath chemicalUnits)
moleking_apply_warnings(berny)

# Collect the frequently reused libraries in one list.
set(MOLEKING_CORE_LIBS myMath chemicalUnits outputProcess berny)

set(MOLEKING_LIB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/MoleKing_Libs")

# Keep intermediate static libs in a dedicated build folder so only the Python
# module lands in site-packages (or, if artifacts are copied, they're grouped).
foreach(_lib IN LISTS MOLEKING_CORE_LIBS)
  set_target_properties(${_lib} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${MOLEKING_LIB_OUTPUT_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${MOLEKING_LIB_OUTPUT_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${MOLEKING_LIB_OUTPUT_DIR}"
  )
endforeach()

# Python extension section.
if(Build_Python)
  # Require Python 3.9+ headers and libraries because the module targets them.
  find_package(Python3 3.9 COMPONENTS Interpreter Development.Module QUIET)
  if(NOT Python3_FOUND)
    message(STATUS "Development.Module not found; attempting fallback with Interpreter only")
    find_package(Python3 3.9 COMPONENTS Interpreter REQUIRED)
  endif()
  # Try to use a system pybind11; otherwise fall back to the vendored copy.
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pybind11/CMakeLists.txt")
      message(STATUS "pybind11 package not found; using the vendored copy in ${CMAKE_CURRENT_SOURCE_DIR}/pybind11")
      add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/pybind11")
    else()
      message(FATAL_ERROR "pybind11 was not located. Install pybind11 or add it as a submodule under ${CMAKE_CURRENT_SOURCE_DIR}/pybind11.")
    endif()
  endif()

  # Build the actual Python module from the pybind11 entry point.
  pybind11_add_module(MoleKing MODULE "${MOLEKING_SRC_DIR}/main.cpp")
  target_link_libraries(MoleKing PRIVATE ${MOLEKING_CORE_LIBS})
  target_include_directories(MoleKing PRIVATE ${MOLEKING_PRIVATE_INCLUDE_DIRS})
  target_compile_definitions(MoleKing PRIVATE Build_Python=1 MOLEKING_VERSION="${PROJECT_VERSION}")
  moleking_apply_warnings(MoleKing)
  set_target_properties(MoleKing PROPERTIES
    OUTPUT_NAME "MoleKing"
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
  )
endif()

# Optional CLI executable section (mirrors the legacy main program).
if(MOLEKING_BUILD_CLI)
  add_executable(MoleKing_cli "${MOLEKING_SRC_DIR}/main2.cpp")
  target_link_libraries(MoleKing_cli PRIVATE ${MOLEKING_CORE_LIBS})
  target_include_directories(MoleKing_cli PRIVATE ${MOLEKING_PRIVATE_INCLUDE_DIRS})
  target_compile_definitions(MoleKing_cli PRIVATE MOLEKING_VERSION="${PROJECT_VERSION}")
  moleking_apply_warnings(MoleKing_cli)
endif()

# Lightweight tests that ensure the CLI reports the expected version.
if(MOLEKING_ENABLE_TESTS)
  enable_testing()
  if(MOLEKING_BUILD_CLI)
    add_test(NAME MoleKingVersion COMMAND MoleKing_cli --version)
    set_tests_properties(MoleKingVersion PROPERTIES PASS_REGULAR_EXPRESSION "MoleKing.*${PROJECT_VERSION}")
  else()
    message(WARNING "MOLEKING_ENABLE_TESTS is ON but MOLEKING_BUILD_CLI is OFF. No executable is available for testing.")
  endif()
endif()
