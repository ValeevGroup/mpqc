set(executable mpqc)

add_executable(
  ${executable}
  mpqc.cc
  mpqcin.cc
  mpqcinit.cc
  parse.cc
  scan.cc
)

target_link_libraries(
  ${executable}
  chemistry
  math
  util
  chemistry
  math
  util
)

target_link_libraries(${executable} ${LIBINT2_LIBRARY})
target_link_libraries(${executable} ${Boost_LIBRARIES})
target_link_libraries(${executable} ${LAPACK_LIBRARIES})
target_link_libraries(${executable} ${BLAS_LIBRARIES})
target_link_libraries(${executable} ${PYTHON_LIBRARIES})

target_link_libraries(${executable} ${LIBRARIES})
target_link_libraries(${executable} ${CMAKE_DL_LIBS})

install(TARGETS ${executable} DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

