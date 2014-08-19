set(matrix_source
  lr_tile.h
  lr_tile.cpp
)

add_library(matrix ${matrix_source})
target_link_libraries(matrix ${TBB_LIBRARIES})

add_executable(matrix_test matrix_test.cpp)
target_link_libraries(matrix_test matrix)
target_link_libraries(matrix_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(matrix_test ${MPI_CXX_LIBRARIES})
