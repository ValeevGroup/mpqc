include(AddLibrary)

# scmat library

add_mpqc_object_library(scmat
  abstract.cc
  blas.cc
  blkiter.cc
  block.cc
  blocked.cc
  blockeddiag.cc
  blockedrect.cc
  blockedsymm.cc
  blockedvect.cc
  c_svd.c
  cmatrix.cc
  dim.cc
  dist.cc
  distdiag.cc
  disthql.cc
  distrect.cc
  distsymm.cc
  distvect.cc
  elemop.cc
  local.cc
  localdiag.cc
  localrect.cc
  localsymm.cc
  localvect.cc
  matrix.cc
  matrix3.cc
  matrix_i.cc
  pdsteqr.F
  repl.cc
  repldiag.cc
  replrect.cc
  replsymm.cc
  replvect.cc
  result.cc
  stateutil.cc
  svd.cc
  util.cc
  vector3.cc
  vector3_i.cc)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE blockedtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(blockedtest blockedtest.cc $<TARGET_OBJECTS:scmat>)
  target_link_libraries(blockedtest
    math
  )
  add_test(blockedtest blockedtest)

  set_property(SOURCE disttest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(disttest disttest.cc $<TARGET_OBJECTS:scmat>)
  target_link_libraries(disttest
    math
  )
  add_test(disttest disttest)

  set_property(SOURCE localtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(localtest localtest.cc $<TARGET_OBJECTS:scmat>)
  target_link_libraries(localtest
    math
  )
  add_test(localtest localtest)

  set_property(SOURCE matrixtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(matrixtest matrixtest.cc $<TARGET_OBJECTS:scmat>)
  target_link_libraries(matrixtest
    math
  )
  add_test(matrixtest matrixtest)

  set_property(SOURCE repltest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(repltest repltest.cc $<TARGET_OBJECTS:scmat>)
  target_link_libraries(repltest
    math
  )
  add_test(repltest repltest)

endif()