include(AddLibrary)

# scf library

add_mpqc_object_library(scf
  clhf.cc
  clscf.cc
  cadf/cadfclhf.cc
  cadf/get_ints.cc
  cadf/decomp.cc
  cadf/coefficients.cc
  cadf/compute_j.cc
  cadf/compute_k.cc
  effh.cc
  fbclhf.cc
  hsoshf.cc
  hsosscf.cc
  osshf.cc
  ossscf.cc
  scf.cc
  scfgradient.cc
  scfops.cc
  scfvector.cc
  svd.cc
  tchf.cc
  tcscf.cc
  tmplinst.cc
  uhf.cc
  uscf.cc
)

# tests

set_property(SOURCE focktest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
add_executable(focktest EXCLUDE_FROM_ALL 
    focktest.cc 
    #$<TARGET_OBJECTS:scf>
)

set_property(SOURCE scftest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
add_executable(scftest EXCLUDE_FROM_ALL 
    scftest.cc 
    #$<TARGET_OBJECTS:scf>
)

target_link_libraries( focktest chemistry math util chemistry math util )
target_link_libraries( scftest chemistry math util chemistry math util )

if (HAVE_LIBINT2)
  target_link_libraries(scftest ${LIBINT2_LIBRARY})
  target_link_libraries(focktest ${LIBINT2_LIBRARY})
else()
  message( FATAL_ERROR "Need libint2" )
endif()

target_link_libraries(scftest ${LAPACK_LIBRARIES})
target_link_libraries(scftest ${BLAS_LIBRARIES})
target_link_libraries(scftest ${LIBRARIES})
target_link_libraries(scftest ${CMAKE_DL_LIBS})

target_link_libraries(focktest ${LAPACK_LIBRARIES})
target_link_libraries(focktest ${BLAS_LIBRARIES})
target_link_libraries(focktest ${LIBRARIES})
target_link_libraries(focktest ${CMAKE_DL_LIBS})
