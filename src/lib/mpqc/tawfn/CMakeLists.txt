set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -std=c++11 -stdlib=libc++")

if (MPQC_UNITTEST)
  add_executable(tawfn_test tawfn_test.cpp )
  target_link_libraries(tawfn_test 
    chemistry
    math
    util
    ${Boost_UNITTEST_LIBRARIES}
  )
  add_test(ta_wavefunction tawfn_test)
endif()
