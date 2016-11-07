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
add_test(build_${executable} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${executable})

add_test(NAME run_${executable}
         COMMAND ${executable})
set_tests_properties(run_${executable}
                     PROPERTIES DEPENDS build_${executable}
                                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
