set(sources
    integral_example.cpp
    )
add_executable(integral_test EXCLUDE_FROM_ALL integral_example.cpp)
target_link_libraries(integral_test ${TileClusterChem_LIBRARIES})
