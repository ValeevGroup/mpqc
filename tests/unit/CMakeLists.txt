enable_testing(true)
add_custom_target(check_unit COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(utests_src
    atom_test.cpp
    clustering_test.cpp
    f12_utility_test.cpp
    keyval_test.cpp
    libint_test.cpp
    orbital_index_test.cpp
    formula_test.cpp
    util_string.cpp
    )

#TODO install xyz file to bin directory

set(executable unit_tests)
add_executable(${executable} EXCLUDE_FROM_ALL test_main.cpp ${utests_src})
target_link_libraries(${executable} ${TA_LIBRARIES} ${Boost_LIBRARIES} ${LIBINT2_LIBRARIES} MPQCbasis MPQCmolecule MPQCexpression MPQCutil_keyval MPQCf12utility)
add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${executable})

add_test(NAME unit/run
         COMMAND ${executable})
set_tests_properties(unit/run
                     PROPERTIES DEPENDS unit/build
                                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/unit
                                ENVIRONMENT MAD_NUM_THREADS=2)
