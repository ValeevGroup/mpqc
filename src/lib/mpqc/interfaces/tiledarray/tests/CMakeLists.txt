# Test of Drew's New Stuff


set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -std=c++11 -stdlib=libc++")

add_executable(molecule_test
   EXCLUDE_FROM_ALL molecule_test.cpp )

target_link_libraries(molecule_test chemistry)
target_link_libraries(molecule_test math)
target_link_libraries(molecule_test util)
target_link_libraries(molecule_test ${LINK_LIBRARIES})


add_executable(basis_test
   EXCLUDE_FROM_ALL basis_test.cpp )

target_link_libraries(basis_test chemistry)
target_link_libraries(basis_test math)
target_link_libraries(basis_test util)
target_link_libraries(basis_test ${LINK_LIBRARIES})


add_executable(integral_test
   EXCLUDE_FROM_ALL integral_test.cpp )

target_link_libraries(integral_test chemistry)
target_link_libraries(integral_test math)
target_link_libraries(integral_test util)
target_link_libraries(integral_test ${LINK_LIBRARIES})


add_executable(fock_build
   EXCLUDE_FROM_ALL fock_build.cpp )

target_link_libraries(fock_build chemistry)
target_link_libraries(fock_build math)
target_link_libraries(fock_build util)
target_link_libraries(fock_build ${LINK_LIBRARIES})


add_executable(fock_build_df
   EXCLUDE_FROM_ALL fock_build_df.cpp )

target_link_libraries(fock_build_df chemistry)
target_link_libraries(fock_build_df math)
target_link_libraries(fock_build_df util)
target_link_libraries(fock_build_df ${LINK_LIBRARIES})


add_executable(dmm_scf
   EXCLUDE_FROM_ALL dmm_scf.cpp )

target_link_libraries(dmm_scf chemistry)
target_link_libraries(dmm_scf math)
target_link_libraries(dmm_scf util)
target_link_libraries(dmm_scf ${LINK_LIBRARIES})