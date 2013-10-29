add_executable(ci.restricted.test EXCLUDE_FROM_ALL test.cpp)
target_link_libraries(ci.restricted.test
  mpqc_test
  ${ARMCI_LIBRARIES}
  ${MPI_LIBRARIES}
)
add_test(ci.restricted ci.restricted.test)
