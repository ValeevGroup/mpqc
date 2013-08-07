set(executable mpqc)

add_library(
  mpqcinit
  SHARED
  mpqcin.cc
  mpqcinit.cc
  parse.cc
  scan.cc  
)

target_link_libraries(
  mpqcinit
  chemistry
  util
)

add_executable(
  ${executable}
  mpqc.cc
)

target_link_libraries(
  ${executable}
  mpqcinit
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

install(TARGETS mpqcinit DESTINATION lib)
