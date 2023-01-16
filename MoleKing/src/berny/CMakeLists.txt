add_library(berny Hessian.cpp)

target_include_directories(berny PRIVATE 
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../myMath"
        "${CMAKE_CURRENT_SOURCE_DIR}/../chemicalUnits"
        )
target_link_libraries(berny
                        PRIVATE myMath
                        PRIVATE chemicalUnits
)

#target_include_directories(berny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")