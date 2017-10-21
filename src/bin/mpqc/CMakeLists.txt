##################### MPQCmpqc library ###########################
set(mpqc_sources
  mpqc_init.cpp
)
set(mpqc_headers
  mpqc_init.h
  mpqc_task.h
)

add_mpqc_library(mpqc mpqc_sources mpqc_headers
                 "MPQCutil_keyval;MPQCutil_mad;MPQCutil_options;MPQCproperties" "mpqc")
target_include_directories(MPQCmpqc PUBLIC ${LIBINT2_INCLUDE_DIRS})
target_link_libraries(MPQCmpqc PUBLIC ${LIBINT2_LIBRARIES})
target_compile_definitions(MPQCmpqc PUBLIC ${LIBINT2_EXTRA_DEFINITIONS})

################## libmpqc super-library #########################

# add the libraries you need and force their linkage in mpqc.cpp
# developers will want to trim this list down to the minimum needed
#set(MPQC_LIBS MPQClcao_f12 MPQClcao_cc MPQClcao_mbpt MPQClcao_scf MPQCmra)
set(MPQC_LIBS MPQClcao_f12 MPQClcao_ci MPQClcao_cc MPQClcao_mbpt MPQClcao_scf)

# create libmpqc, an interface to MPQC_LIBS
# linking libmpqc in will link in ALL MPQC functionality
set(deplibs "MPQCmpqc" ${MPQC_LIBS})
add_library(libmpqc INTERFACE)
target_link_libraries(libmpqc INTERFACE ${deplibs})
install(TARGETS libmpqc EXPORT mpqc
        COMPONENT libmpqc
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")

add_custom_target(install-libmpqc)
add_dependencies(install-libmpqc install-mpqc)

################## MPQCmain useful to develop plugins #########################
set(main_source mpqc.cpp)
add_mpqc_library(main  main_source ""
        "MPQCmpqc" "mpqc")
target_link_libraries(MPQCmain PUBLIC ${MPQC_LIBS})
# should make MPQCmain a MODULE ?? it's missing ForceLink symbols, but that's up to the user to satisfy
if(APPLE)
#    target_link_libraries(MPQCmain PUBLIC "-undefined dynamic_lookup")
endif(APPLE)

##################### MPQC executable ###########################
add_executable(mpqc mpqc.cpp)
INCLUDE_DIRECTORIES(${MPQC_SOURCE_DIR}/src/bin)
target_link_libraries(mpqc
                      libmpqc ${MADNESS_DISABLEPIE_LINKER_FLAG})
# note above an ad hoc treatment of PIE until CMake allows to do it properly ... specifically, need CMake to:
# 1. force non-PIE for compilers that default to PIE
# 2. allow clean specification of PIE/non-PIE properties for library targets

install(TARGETS mpqc RUNTIME DESTINATION bin)
#cotire(mpqc)

add_executable(mpqc_test mpqc_test.cpp)
set_property(SOURCE mpqc_test.cpp PROPERTY COMPILE_DEFINITIONS
             SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mpqc_test libmpqc ${MADNESS_DISABLEPIE_LINKER_FLAG})
add_test(mpqc_test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target mpqc_test)

# multiworld test is broken now, only run 1 MPI rank for now
add_test(NAME mpqc_test/run
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} mpqc_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(mpqc_test/run
                     PROPERTIES DEPENDS mpqc_test/build
                                ENVIRONMENT MAD_NUM_THREADS=2)

