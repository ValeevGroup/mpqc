add_custom_target(examples)

add_executable(mpqc_info mpqc_info.cpp)
target_link_libraries(mpqc_info MPQCbasis MPQCmolecule MPQCutil_keyval MPQCexpression MPQCmath_ta MPQCutil_misc ${Boost_LIBRARIES})

add_executable(phf EXCLUDE_FROM_ALL phf.cpp)
target_link_libraries(phf MPQCmolecule MPQCkeyval MPQCmath_ta MPQCbasis MPQCexpression MPQCintegrals)

