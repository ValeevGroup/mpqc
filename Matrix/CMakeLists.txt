set(SOURCES
    full_rank_tile.h
    lin_algebra.h
    low_rank_tile.h
    lr_tile.h
    tile_concept.h
)

add_executable(matrix_test matrix_test.cpp)
target_link_libraries(matrix_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(matrix_test ${LAPACK_LIBRARIES})
target_link_libraries(matrix_test ${MPI_CXX_LIBRARIES})

add_executable(lrtile_reduced_only_test lrtile_reduced_only_test.cpp)
target_link_libraries(lrtile_reduced_only_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(lrtile_reduced_only_test ${LAPACK_LIBRARIES})
target_link_libraries(lrtile_reduced_only_test ${MPI_CXX_LIBRARIES})

add_executable(lrtile_full_only_test lrtile_full_only_test.cpp)
target_link_libraries(lrtile_full_only_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(lrtile_full_only_test ${LAPACK_LIBRARIES})
target_link_libraries(lrtile_full_only_test ${MPI_CXX_LIBRARIES})

add_executable(lr_reduced_full_test lr_reduced_full_test.cpp)
target_link_libraries(lr_reduced_full_test ${TILEDARRAY_LIBRARIES})
target_link_libraries(lr_reduced_full_test ${LAPACK_LIBRARIES})
target_link_libraries(lr_reduced_full_test ${MPI_CXX_LIBRARIES})

add_subdirectory(eigen_tests)
