set(executable mpqc)

add_library(
  mpqcinit
  mpqcin.cc
  mpqcinit.cc
  parse.cc
  scan.cc  
)

add_executable(
  ${executable}
  mpqc.cc
)

target_link_libraries(
  ${executable}
  mpqcinit
  chemistry
  math
  util
  chemistry
  math
  util
)

if (HAVE_LIBINT2)
  target_link_libraries(${executable} ${LIBINT2_LIBRARY})
endif()

target_link_libraries(${executable} ${Boost_LIBRARIES})
target_link_libraries(${executable} ${LAPACK_LIBRARIES})
target_link_libraries(${executable} ${BLAS_LIBRARIES})
target_link_libraries(${executable} ${PYTHON_LIBRARIES})

target_link_libraries(${executable} ${LIBRARIES})
target_link_libraries(${executable} ${CMAKE_DL_LIBS})

install(TARGETS ${executable} DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

set(headers mpqcinit.h mpqcin.h)
foreach (_file ${headers})
  #message(STATUS "gonna install ${_file} in include/mpqc" )
  install(FILES ${_file} DESTINATION include/mpqc)
endforeach()

install(TARGETS mpqcinit ARCHIVE DESTINATION lib)
