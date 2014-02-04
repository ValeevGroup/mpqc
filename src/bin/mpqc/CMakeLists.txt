set(executable mpqc)

add_library(
  mpqcmain
  mpqc.cc
)

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
  mpqcmain
  mpqcinit
)

target_link_libraries(
  mpqcinit
  chemistry
  util
)

# N.B. mpqcinit alone has almost no dependencies
# library dependencies should be tied to chemistry and util
target_link_libraries(
  ${executable}
  mpqcinit
)

set_target_properties(${executable} PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY")
cotire(${executable})

# message("** LINK_LIBRARIES=${LINK_LIBRARIES}")

install(TARGETS ${executable} RUNTIME DESTINATION bin)

set(headers mpqcinit.h mpqcin.h)
foreach (_file ${headers})
  #message(STATUS "gonna install ${_file} in include/mpqc" )
  install(FILES ${_file} DESTINATION include/mpqc)
endforeach()

install(TARGETS mpqcinit DESTINATION lib)
install(TARGETS mpqcmain DESTINATION lib)
