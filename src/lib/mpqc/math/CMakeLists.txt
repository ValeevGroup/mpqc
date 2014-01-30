if(MPQC_UNITTEST)
  add_executable(blas_test blas_test.cpp)

  target_link_libraries(blas_test
    util
    ${LAPACK_LIBRARIES}
    ${Boost_UNITTEST_LIBRARIES}
  )

  add_test(blas blas_test)
endif()
