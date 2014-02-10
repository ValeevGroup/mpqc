include(AddLibrary)

# wfn library

# collect the .cc and .cpp files to be compiled.
set(ccfiles 
  accum.cc
  density.cc
  densval.cc
  eht.cc
  esp.cc
  femo.cc
  hcorewfn.cc
  nao.cc
  natorbval.cc
  obwfn.cc
  orbital.cc
  orbitalspace.cc
  orbitalspace_utils.cc
  rdm.cc
  spin.cc
  wfn.cc
)

#add new features if present
IF(MPQC_NEW_FEATURES)
  LIST(APPEND ccfiles tawfn.cpp)
ENDIF(MPQC_NEW_FEATURES)

add_mpqc_object_library(wfn ${ccfiles})

# tests

#set_property(SOURCE wfntest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(wfntest EXCLUDE_FROM_ALL wfntest.cc $<TARGET_OBJECTS:wfn>)
#add_test(wfntest wfntest)

