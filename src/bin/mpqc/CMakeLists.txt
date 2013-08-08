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

# Ed, the previous linking configuration broke badly
# afaik, only one target_link_libraries per target is allowed
target_link_libraries(
  ${executable}
  mpqcinit
  chemistry
  util
  ${LINK_LIBRARIES})

#message("** LINK_LIBRARIES=${LINK_LIBRARIES}")

install(TARGETS ${executable} RUNTIME DESTINATION bin)

set(headers mpqcinit.h mpqcin.h)
foreach (_file ${headers})
  #message(STATUS "gonna install ${_file} in include/mpqc" )
  install(FILES ${_file} DESTINATION include/mpqc)
endforeach()

install(TARGETS mpqcinit DESTINATION lib)
