add_executable(tiledbasis_test tiledbasisset_test.cpp)



target_link_libraries(tiledbasis_test
    chemistry
    ${Boost_LIBRARIES}
)

add_test(trange1 trange1_test)