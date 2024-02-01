#add_subdirectory(../myMath)
add_library(chemicalUnits 
        Molecule.hpp
        Molecule.cpp
        AtomicScale.hpp
        AtomicScale.cpp
        PeriodicTable.hpp
        PeriodicTable.cpp
        SupraMolecule.hpp
        SupraMolecule.cpp
        OPLSff.hpp
        OPLSff.cpp
        pov.cpp
        pov.hpp) 
target_include_directories(chemicalUnits PRIVATE 
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../myMath"
        )
target_link_libraries(chemicalUnits
                        PRIVATE myMath
        )
#target_include_directories(chemicalUnits PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")