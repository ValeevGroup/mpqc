add_executable(blas_test  EXCLUDE_FROM_ALL blas_test.cpp)

target_link_libraries(blas_test
  ${LAPACK_LIBRARIES}
  ${Boost_TEST_LIBRARY}
  ${Boost_LIBRARIES}
)

add_test(blas blas_test)
