include(AddLibrary)

# libint2 library

add_mpqc_object_library(libint2
  bounds.cc
  comp_eri.cc
  comp_g12.cc
  comp_g12dkh.cc
  comp_g12nc.cc
  #comp_grt.cc
  edipole.cc
  equadrupole.cc
  eri.cc
  g12.cc
  g12dkh.cc
  g12nc.cc
  #grt.cc
  gto.cc
  hcore.cc
  int1e.cc
  int2e.cc
  kinetic.cc
  libint2.cc
  nuclear.cc
  obintlibint2.cc
  obosrr.cc
  overlap.cc
  p4.cc
  permute2e.cc
  primpairs.cc
  shellpairs.cc
  storage.cc
  tbintlibint2.cc
  tform.cc
)

# tests

set_property(SOURCE libint2test.cc PROPERTY COMPILE_DEFINITIONS
   SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(libint2test EXCLUDE_FROM_ALL
 libint2test.cc
)

target_link_libraries(
  libint2test
  chemistry
)

target_link_libraries(libint2test ${LIBINT2_LIBRARY})
target_link_libraries(libint2test ${LAPACK_LIBRARIES})
target_link_libraries(libint2test ${BLAS_LIBRARIES})
target_link_libraries(libint2test ${LIBRARIES})
target_link_libraries(libint2test ${CMAKE_DL_LIBS})

