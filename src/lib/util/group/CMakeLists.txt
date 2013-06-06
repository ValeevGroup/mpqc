# group library

set(sources
  actmsg.cc
#  file.cc
#  fileproc.cc
  globcnt.cc
  hcube.cc
  memamsg.cc
#  memarmci.cc
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

add_library(group OBJECT ${sources})

# tests

set_property(SOURCE actmsgtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(actmsgtest EXCLUDE_FROM_ALL actmsgtest.cc $<TARGET_OBJECTS:group>)
add_test(actmsgtest actmsgtest)

set_property(SOURCE memtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(memtest EXCLUDE_FROM_ALL memtest.cc $<TARGET_OBJECTS:group>)
add_test(memtest memtest)

set_property(SOURCE messtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(messtest EXCLUDE_FROM_ALL messtest.cc $<TARGET_OBJECTS:group>)
add_test(messtest messtest)

set_property(SOURCE pooltest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(pooltest EXCLUDE_FROM_ALL pooltest.cc $<TARGET_OBJECTS:group>)
add_test(pooltest pooltest)

set_property(SOURCE prttest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(prttest EXCLUDE_FROM_ALL prttest.cc $<TARGET_OBJECTS:group>)
add_test(prttest prttest)

set_property(SOURCE rnglocktest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(rnglocktest EXCLUDE_FROM_ALL rnglocktest.cc $<TARGET_OBJECTS:group>)
add_test(rnglocktest rnglocktest)

set_property(SOURCE thrtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(thrtest EXCLUDE_FROM_ALL thrtest.cc $<TARGET_OBJECTS:group>)
add_test(thrtest thrtest)

