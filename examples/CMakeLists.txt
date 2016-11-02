add_custom_target(examples)

add_executable(mpqc_test EXCLUDE_FROM_ALL mpqc_test.cpp)
target_link_libraries(mpqc_test ${MPQC_LIBRARIES} MPQCmisc)

add_executable(mpqc_info mpqc_info.cpp)
target_link_libraries(mpqc_info MPQCbasis MPQCmolecule MPQCkeyval MPQCexpression MPQCmath_ta MPQCmisc ${Boost_LIBRARIES})

