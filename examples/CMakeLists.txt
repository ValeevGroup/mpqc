add_custom_target(examples)

add_executable(mpqc_info mpqc_info.cpp)
target_link_libraries(mpqc_info MPQCbasis MPQCmolecule MPQCkeyval MPQCexpression MPQCmath_ta MPQCmisc ${Boost_LIBRARIES})

add_executable(phf EXCLUDE_FROM_ALL phf.cpp)
target_link_libraries(phf ${MPQC_LIBRARIES})

