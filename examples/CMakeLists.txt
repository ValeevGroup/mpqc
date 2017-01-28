add_custom_target(examples)

add_executable(mpqc_info mpqc_info.cpp)
target_link_libraries(mpqc_info MPQClcao_basis MPQCmolecule MPQCutil_keyval MPQClcao_expression MPQCmath_ta MPQCutil_misc MPQCutil_mad ${Boost_LIBRARIES})

