add_executable(tile_algebra_examples tile_algebra.cpp)
target_link_libraries(tile_algebra_examples
    ${TILEDARRAY_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)
