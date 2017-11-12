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
# libmpqc is an interface to the *entire* top-level MPQC functionality
# this target is useful to create a standalone executable that augments the MPQC functionality
# N.B. To force linking in of these libs into the executqble use ForceLink mechanism (see the top of mpqc.cpp)
# N.B. due to transitive dependencies this could be shortened to: set (mpqc_features "lcao_f12;lcao_ci")
#      *but* this would break introspection of MPQC capabilities later (e.g. in validation suite)
set(mpqc_features "lcao_f12;lcao_cc;lcao_mbpt;lcao_ci;lcao_scf")
# N.B. developers will want to trim this list down to the minimum needed, e.g. down to "lcao_scf" to test HF only
#set (mpqc_features "lcao_scf" )
set(MPQC_FEATURES ${mpqc_features} CACHE STRING
    "The list of features supported by the MPQC executable" FORCE)
# transform MPQC_FEATURES to the list of library targets and associated preprocessor definitions
set(MPQC_FEATURE_LIBS )
set(MPQC_FEATURE_DEFINITIONS )
foreach(feature ${MPQC_FEATURES})
     list(APPEND MPQC_FEATURE_LIBS "MPQC${feature}")
     list(APPEND MPQC_FEATURE_DEFINITIONS "-DMPQC_INCLUDES_${feature}=1")
endforeach()

# create libmpqc, an interface to MPQC_FEATURE_LIBS
# linking libmpqc in will link in ALL MPQC functionality
set(deplibs "MPQCmpqc" ${MPQC_FEATURE_LIBS})
add_library(libmpqc INTERFACE)
target_link_libraries(libmpqc INTERFACE ${deplibs})
install(TARGETS libmpqc EXPORT mpqc
        COMPONENT libmpqc
        LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}")

add_custom_target(install-libmpqc)
add_dependencies(install-libmpqc install-mpqc)

################## MPQCmain useful to develop plugins #########################
# MPQCmain (see below) also depends on MPQC_FEATURE_LIBS functionality ... can't use libmpqc since it's just an INTERFACE
set(main_source mpqc.cpp)
add_mpqc_library(main  main_source ""
        "MPQCmpqc" "mpqc")
target_link_libraries(MPQCmain PUBLIC ${MPQC_FEATURE_LIBS})
target_compile_definitions(MPQCmain PUBLIC ${MPQC_FEATURE_DEFINITIONS})

##################### MPQC executable ###########################
add_executable(mpqc mpqc.cpp)
INCLUDE_DIRECTORIES(${MPQC_SOURCE_DIR}/src/bin)
target_link_libraries(mpqc
                      libmpqc ${MADNESS_DISABLEPIE_LINKER_FLAG})
target_compile_definitions(mpqc PUBLIC ${MPQC_FEATURE_DEFINITIONS})
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

