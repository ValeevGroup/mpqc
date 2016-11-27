##################### MPQCinit library ###########################
set(init_sources
  mpqc_init.cpp
)
set(init_headers
  mpqc_init.h
)

add_mpqc_library(init init_sources init_headers
                 "MPQCutil_keyval;MPQCutil_mad;MPQCutil_options" "mpqc")

################## libmpqc super-library #########################

# add the libraries you need and force their linkage in mpqc.cpp
# developers will want to trim this list down to the minimum needed
set(MPQC_LIBS MPQCf12 MPQCcc MPQCmbpt MPQCscf)

# create a dummy library MPQCmpqc to be immediately aliased to libmpqc
# linking libmpqc in will link in ALL MPQC functionality
set(libmpqc_sources libmpqc.cpp)
add_mpqc_library(mpqc libmpqc_sources ""
                 ${MPQC_LIBS} "mpqc")
add_library(libmpqc INTERFACE)
target_link_libraries(libmpqc INTERFACE MPQCmpqc)
install(TARGETS libmpqc EXPORT mpqc
        COMPONENT libmpqc
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")

add_custom_target(install-libmpqc)
add_dependencies(install-libmpqc install-mpqc)

##################### MPQC executable ###########################
add_executable(mpqc mpqc.cpp)
                      
target_link_libraries(mpqc
                      MPQCmpqc MPQCinit ${MPQC_LIBS})
install(TARGETS mpqc RUNTIME DESTINATION bin)
#cotire(mpqc)

add_executable(mpqc_test mpqc_test.cpp)
set_property(SOURCE mpqc_test.cpp PROPERTY COMPILE_DEFINITIONS
             SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mpqc_test MPQCinit MPQCscf MPQCutil_options MPQCutil_misc)
add_test(mpqc_test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target mpqc_test)

# multiworld test is broken now, only run 1 MPI rank for now
add_test(NAME mpqc_test/run
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} mpqc_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(mpqc_test/run
                     PROPERTIES DEPENDS mpqc_test/build
                                ENVIRONMENT MAD_NUM_THREADS=2)
                      