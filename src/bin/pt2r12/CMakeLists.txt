
if (GAMESS)

  set(CMAKE_CXX_FLAGS "-DPT2R12GAMESS ${CMAKE_CXX_FLAGS}")
  
endif()

set(executable pt2r12)
add_executable(
    ${executable}
    pt2r12.cc
)

target_link_libraries(
    ${executable}
    mpqcinit
    chemistry
    math
    util
    #chemistry
    #math
    #util
    extern_moinfo
    ${LINK_LIBRARIES}
)

install(TARGETS ${executable} DESTINATION bin)

  #set(headers extern_pt2r12.h moinfo.h molden_moinfo.h)
  #foreach (_file ${headers})
  #  install(FILES ${_file} DESTINATION include/mpqc/pt2r12)
  #endforeach()
