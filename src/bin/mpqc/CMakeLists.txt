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

target_link_libraries(
  ${executable}
  ${LINK_LIBRARIES})

install(TARGETS ${executable} DESTINATION bin)

set(headers mpqcinit.h mpqcin.h)
foreach (_file ${headers})
  #message(STATUS "gonna install ${_file} in include/mpqc" )
  install(FILES ${_file} DESTINATION include/mpqc)
endforeach()

install(TARGETS mpqcinit ARCHIVE DESTINATION lib)
