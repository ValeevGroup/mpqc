# Get GTest
add_subdirectory(gtest-1.7.0)
include_directories(src ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(test_sources 
  full_rank_tile_test.cpp
  low_rank_tile_test.cpp
  basis_set_test.cpp
  tile_algebra_test.cpp
  tile_variant_test.cpp
  tile_ops_test.cpp
  tile_pimpl_test.cpp
  tile_mutations_test.cpp
  all_same_test.cpp
)

FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR} ano_basis_G94.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(unit_tests ${test_sources})
target_link_libraries(unit_tests
  gtest
  gtest_main
  ${TileClusterChem_LIBRARIES}
)

