if(DOES_NOT_EXIST_THINGY)
  add_executable(trange1_test EXCLUDE_FROM_ALL trange1_test.cpp)
  target_link_libraries(trange1_test
    chemistry
    ${Boost_TEST_LIBRARY}
    ${Boost_LIBRARIES}
  )
  add_test(trange1 trange1_test)
endif()