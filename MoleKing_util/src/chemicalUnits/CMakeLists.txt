add_library(chemicalUnits 
        Molecule.hpp
        Molecule.cpp
        AtomicScale.hpp
        AtomicScale.cpp
        PeriodicTable.hpp
        PeriodicTable.cpp
        SupraMolecule.hpp
        SupraMolecule.cpp) 

target_include_directories(chemicalUnits PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")