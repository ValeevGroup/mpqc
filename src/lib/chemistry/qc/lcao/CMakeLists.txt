include(AddLibrary)

# lcao library

set(ccfiles
  clhfcontrib.cc
  compute_iRjS.cc
  compute_ijxy.cc
  compute_ikjy.cc
  compute_ixjy.cc
  df.cc
  df_runtime.cc
  fockbuild.cc
  fockbuild_runtime.cc
  fockbuilder.cc
  fockdist.cc
  hsoshfcontrib.cc
  moints_runtime.cc
  soad.cc
  tbint_runtime.cc
  transform_123inds.cc
  transform_12inds.cc
  transform_13inds.cc
  transform_factory.cc
  transform_iRjS.cc
  transform_ijR.cc
  transform_ijxy.cc
  transform_ikjy.cc
  transform_ixjy.cc
  transform_ixjy_df.cc
  transform_tbint.cc
  utils.cc
  wfnworld.cc
)

add_mpqc_object_library(lcao ${ccfiles})

# tests

set_property(SOURCE dftest.cc PROPERTY COMPILE_DEFINITIONS
   SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(dftest EXCLUDE_FROM_ALL
 dftest.cc
)

target_link_libraries(
  dftest
  chemistry
)

if (HAVE_LIBINT2)
  target_link_libraries(dftest ${LIBINT2_LIBRARY})
endif()

target_link_libraries(dftest ${LAPACK_LIBRARIES})
target_link_libraries(dftest ${BLAS_LIBRARIES})

target_link_libraries(dftest ${LIBRARIES})
target_link_libraries(dftest ${CMAKE_DL_LIBS})

