
if (GAMESS)
  set(CMAKE_CXX_FLAGS "-DMOINFOGAMESS ${CMAKE_CXX_FLAGS}")
endif(GAMESS)

add_library(extern_moinfo
  moinfo.cc
  #molden_moinfo.cc
)
#target_link_libraries(extern_moinfo chemistry)
install(TARGETS extern_moinfo DESTINATION lib)