include(AddLibrary)

# nbody library

SET(ccfiles
  moinfo.cc
)

if (GAMESS)
  set(CMAKE_CXX_FLAGS "-DMOINFOGAMESS ${CMAKE_CXX_FLAGS}")
endif(GAMESS)


add_mpqc_object_library(extern ${ccfiles})
