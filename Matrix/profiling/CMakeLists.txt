add_executable(low_rank_gemm_profile low_rank_gemm.cpp)
target_link_libraries(low_rank_gemm_profile
    ${TILEDARRAY_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

