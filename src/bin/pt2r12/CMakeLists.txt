
if (GAMESS)

  add_library(
    pt2r12gamess
    pt2r12.cc
    moinfo.cc
    extern_pt2r12.cc 
  )

  set(CMAKE_CPP_FLAGS "-DPT2R12GAMESS -DSKIP_MAIN ${CMAKE_CPP_FLAGS}")
  install(TARGETS pt2r12gamess ARCHIVE DESTINATION lib)

else()

  set(executable pt2r12)
  add_executable(
    ${executable}
    pt2r12.cc
    moinfo.cc
    extern_pt2r12.cc
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
    ${LINK_LIBRARIES}
  )

  install(TARGETS ${executable} DESTINATION bin)

  #set(headers extern_pt2r12.h moinfo.h molden_moinfo.h)
  #foreach (_file ${headers})
  #  install(FILES ${_file} DESTINATION include/mpqc/pt2r12)
  #endforeach()

endif()
