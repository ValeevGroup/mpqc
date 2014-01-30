include(AddLibrary)

# lmp2 library

add_mpqc_object_library(lmp2
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
)

# tests
  
if (MPQC_UNITTEST)
  set_property(SOURCE avlmmaptest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(avlmmaptest EXCLUDE_FROM_ALL avlmmaptest.cc $<TARGET_OBJECTS:lmp2>)
  target_link_libraries(avlmmaptest 
    chemistry
  )
  # this test takes too long
  #add_test(avlmmaptest avlmmaptest)

  set_property(SOURCE smatest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(smatest smatest.cc $<TARGET_OBJECTS:lmp2>)
  target_link_libraries(smatest
    chemistry
  )
  add_test(smatest smatest)
  
endif()
