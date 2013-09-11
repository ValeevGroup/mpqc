add_executable(tiledbasis_test EXCLUDE_FROM_ALL tiledbasisset_test.cpp)

target_link_libraries(tiledbasis_test
    chemistry
    ${Boost_TEST_LIBRARY}
    ${Boost_LIBRARIES}
)
