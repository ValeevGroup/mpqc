add_library(
  MPQCinit
  mpqc_init.cpp
  mpqc_init.h
)

target_link_libraries(
  MPQCinit
  MPQCutil_options
  MPQCutil_keyval
  MPQCutil_misc
)

add_executable(mpqc EXCLUDE_FROM_ALL mpqc.cpp)
# add the libraries you need and force their linkage in mpqc.cpp
# developers will want to trim this list down to the minimum needed
set(MPQC_LIBS MPQCf12 MPQCcc MPQCmbpt MPQCscf MPQCphf)
                      
target_link_libraries(mpqc
                      ${MPQC_LIBS}
                      MPQCutil_options MPQCutil_misc
                      MPQCinit)

add_executable(mpqc_test mpqc_test.cpp)
set_property(SOURCE mpqc_test.cpp PROPERTY COMPILE_DEFINITIONS
             SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mpqc_test MPQCscf MPQCutil_options MPQCutil_misc
                      MPQCinit)
add_test(mpqc_test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target mpqc_test)

# multiworld test is broken now, only run 1 MPI rank for now
add_test(NAME mpqc_test/run
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} mpqc_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(mpqc_test/run
                     PROPERTIES DEPENDS mpqc_test/build
                                ENVIRONMENT MAD_NUM_THREADS=2)

