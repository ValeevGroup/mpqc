# Get GTest
add_subdirectory(gtest-1.7.0)
include_directories(src ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(test_sources 
  full_rank_tile_test.cpp
)

add_executable(unit_tests full_rank_tile_test.cpp)
target_link_libraries(unit_tests
  gtest
  gtest_main
  ${HeirChem_LIBRARIES}
)

