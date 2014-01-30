include(AddLibrary)

# optimize library

add_mpqc_object_library(optimize
  conv.cc
  dfp.cc
  diis.cc
  efc.cc
  function.cc
  gaussianfit.cc
  gdiis.cc
  lbfgs.f
  mcsearch.cc
  newton.cc
  opt.cc
  powell.cc
  qnewton.cc
  scextrap.cc
  scextrapmat.cc
  steep.cc
  transform.cc
  update.cc
  levmar/lm.c
  levmar/Axb.c
  levmar/lmbc.c
  levmar/misc.c
)

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE opttest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(opttest opttest.cc $<TARGET_OBJECTS:optimize>)
  target_link_libraries(opttest
    math
  )
  add_test(opttest opttest)

  set_property(SOURCE scextest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(scextest scextest.cc $<TARGET_OBJECTS:optimize>)
  target_link_libraries(scextest
    math
  )
  add_test(scextest scextest)
endif()

