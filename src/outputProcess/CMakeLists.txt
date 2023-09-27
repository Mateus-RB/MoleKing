add_library(outputProcess 
                        G16Process.cpp
                        G16Process.hpp
                        Psi4Process.cpp
                        Psi4Process.hpp)

target_include_directories( outputProcess PRIVATE 
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../myMath"
        "${CMAKE_CURRENT_SOURCE_DIR}/../chemicalUnits"
        )
target_link_libraries(outputProcess
                        PRIVATE myMath
                        PRIVATE chemicalUnits
                    )

#target_include_directories(outputProcess PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")