include(AddLibrary)

# basis library

add_mpqc_object_library(basis
  aotoso.cc
  cart.cc
  cartiter.cc
  dercent.cc
  distsh.cc
  distshpair.cc
  extent.cc
  files.cc
  fjt.cc
  gaussbas.cc
  gaussbaval.cc
  gaussshell.cc
  gaussshval.cc
  gpetite.cc
  integral.cc
  intparams.cc
  inttraits.cc
  lselect.cc
  obint.cc
  operator.cc
  orthog.cc
  petite.cc
  shellrot.cc
  sobasis.cc
  sointegral.cc
  split.cc
  symmint.cc
  tbint.cc
  transform.cc
  uncontract.cc
  union.cc
)

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE btest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(btest btest.cc $<TARGET_OBJECTS:basis>)
  add_test(btest btest)

  target_link_libraries( btest
    chemistry
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LIBRARIES}
    ${CMAKE_DL_LIBS}
  )
endif()
