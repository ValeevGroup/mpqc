add_executable(low_rank_gemm_profile low_rank_gemm.cpp)
target_link_libraries(low_rank_gemm_profile
    ${TILEDARRAY_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

add_executable(low_rank_subt_to_profile low_rank_subt_to.cpp)
target_link_libraries(low_rank_subt_to_profile
    ${TILEDARRAY_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

add_executable(full_rank_gemm_profile full_rank_gemm.cpp)
target_link_libraries(full_rank_gemm_profile
    ${TILEDARRAY_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

