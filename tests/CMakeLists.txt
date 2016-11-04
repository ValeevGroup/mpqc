set(utests
    atom_test.cpp
    clustering_test.cpp
    f12_utility_test.cpp
    keyval_test.cpp
    orbital_index_test.cpp
    formula_test.cpp
    util_string.cpp
    )

#TODO install xyz file to bin directory

add_executable(tests EXCLUDE_FROM_ALL test_main.cpp ${utests})
target_link_libraries(tests ${TA_LIBRARIES} ${Boost_LIBRARIES} MPQCbasis MPQCmolecule MPQCexpression MPQCutil_keyval MPQCf12utility)
add_test(NAME TestAll COMMAND ${CMAKE_BINARY_DIR}/tests/tests)
