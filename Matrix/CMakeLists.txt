add_executable(matrix_test matrix_test.cpp)
target_link_libraries(matrix_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(matrix_test ${LAPACK_LIBRARIES})
target_link_libraries(matrix_test ${MPI_CXX_LIBRARIES})

add_executable(LRTile_test LRTile_test.cpp)
target_link_libraries(LRTile_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(LRTile_test ${LAPACK_LIBRARIES})
target_link_libraries(LRTile_test ${MPI_CXX_LIBRARIES})

add_subdirectory(eigen_tests)
