enable_testing(true)
add_custom_target(check_unit COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(utests_src
    atom_test.cpp
    bug_test.cpp
    clustering_test.cpp
        davidson_diag_test.cpp
    exception_test.cpp
    f12_utility_test.cpp
    formio_test.cpp
    formula_test.cpp
        gram_schmidt_test.cpp
    keyval_test.cpp
    libint_test.cpp
    molecule_test.cpp
    orbital_index_test.cpp
    units_test.cpp
    util_string.cpp
    )

#TODO install xyz file to bin directory

set(executable unit_tests)
add_executable(${executable} EXCLUDE_FROM_ALL test_main.cpp ${utests_src})
target_link_libraries(${executable} MPQCmpqc MPQClcao_integrals)
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${executable})
set_property(SOURCE bug_test.cpp PROPERTY COMPILE_DEFINITIONS
             SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME unit/run-np-1
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${executable}> ${MPIEXEC_POSTFLAGS})
set_tests_properties(unit/run-np-1
                     PROPERTIES DEPENDS unit/build
                                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit
                                ENVIRONMENT MAD_NUM_THREADS=2)
add_test(NAME unit/run-np-2
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${executable}> ${MPIEXEC_POSTFLAGS})
set_tests_properties(unit/run-np-2
                     PROPERTIES DEPENDS unit/build
                                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit
                                ENVIRONMENT MAD_NUM_THREADS=2)
