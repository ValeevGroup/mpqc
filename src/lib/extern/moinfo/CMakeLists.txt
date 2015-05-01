
add_library(extern_moinfo
  moinfo.cc
  #molden_moinfo.cc
)
#target_link_libraries(extern_moinfo chemistry)
install(TARGETS extern_moinfo DESTINATION lib)