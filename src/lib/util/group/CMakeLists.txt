include(AddLibrary)

# group library

set(sources
  actmsg.cc
#  file.cc
#  fileproc.cc
  globcnt.cc
  hcube.cc
  memamsg.cc
  memiter.cc
  memmsg.cc
  memory.cc
  memproc.cc
  memrdma.cc
  memregion.cc
  memshm.cc
  message.cc
  messimpl.cc
  messproc.cc
  mstate.cc
  pool.cc
  pregtime.cc
  reduce.cc
  rnglock.cc
  thpthd.cc
#  thpuma.cc
  thread.cc
  topology.cc)

if(HAVE_MPI)
  list(APPEND sources messmpi.cc memmtmpi.cc)
endif()

# # ARMCIMemGroup broken, doesnt compile
# if (HAVE_ARMCI)
#   list(APPEND sources memarmci.cc)
# endif()

add_mpqc_object_library(group ${sources})

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE actmsgtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(actmsgtest actmsgtest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(actmsgtest
    util
  )
  add_test(actmsgtest actmsgtest)

  set_property(SOURCE memtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(memtest memtest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(memtest
    util
  )
  add_test(memtest memtest)

  set_property(SOURCE messtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(messtest messtest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(messtest
    util
  )
  #add_test(messtest messtest)

  set_property(SOURCE pooltest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(pooltest pooltest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(pooltest
    util
  )
  add_test(pooltest pooltest)

  set_property(SOURCE prttest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(prttest prttest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(prttest
    util
  )
  add_test(prttest prttest)

  set_property(SOURCE rnglocktest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(rnglocktest rnglocktest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(rnglocktest
    util
  )
  add_test(rnglocktest rnglocktest)

  set_property(SOURCE thrtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(thrtest thrtest.cc $<TARGET_OBJECTS:group>)
  target_link_libraries(thrtest
    util
  )
  add_test(NAME thrtest COMMAND thrtest -threadgrp "<PthreadThreadGrp>:(num_threads=2)")

endif()