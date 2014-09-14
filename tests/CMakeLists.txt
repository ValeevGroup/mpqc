# Get GTest
add_subdirectory(gtest-1.7.0)
include_directories(src ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(test_sources 
  full_rank_tile_test.cpp
  low_rank_tile_test.cpp
  tile_algebra_test.cpp
  tileable_test.cpp
)

add_executable(unit_tests ${test_sources})
target_link_libraries(unit_tests
  gtest
  gtest_main
  ${TileClusterChem_LIBRARIES}
)

