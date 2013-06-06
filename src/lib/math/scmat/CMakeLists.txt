# scmat library

add_library(scmat OBJECT
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
  pdsteqr.f
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

set_property(SOURCE blockedtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(blockedtest EXCLUDE_FROM_ALL blockedtest.cc $<TARGET_OBJECTS:scmat>)
add_test(blockedtest blockedtest)

set_property(SOURCE disttest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(disttest EXCLUDE_FROM_ALL disttest.cc $<TARGET_OBJECTS:scmat>)
add_test(disttest disttest)

set_property(SOURCE localtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(localtest EXCLUDE_FROM_ALL localtest.cc $<TARGET_OBJECTS:scmat>)
add_test(localtest localtest)

set_property(SOURCE matrixtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(matrixtest EXCLUDE_FROM_ALL matrixtest.cc $<TARGET_OBJECTS:scmat>)
add_test(matrixtest matrixtest)

set_property(SOURCE repltest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(repltest EXCLUDE_FROM_ALL repltest.cc $<TARGET_OBJECTS:scmat>)
add_test(repltest repltest)

