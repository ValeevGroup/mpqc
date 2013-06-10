# lmp2 library

add_library(lmp2 OBJECT
  collective.cc
  dgemminfo.cc
  domain.cc
  extrap.cc
  lcorr.cc
  lmp2.cc
  pack2e.cc
  pop_local.cc
  sma.cc
  smautil.cc
  testsma.cc
  timesma.cc
)

# tests

#set_property(SOURCE avlmmaptest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(avlmmaptest EXCLUDE_FROM_ALL avlmmaptest.cc $<TARGET_OBJECTS:lmp2>)
#add_test(avlmmaptest avlmmaptest)

