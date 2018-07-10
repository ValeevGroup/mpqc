add_custom_target(examples)

add_executable(mpqc_info mpqc_info.cpp)
target_link_libraries(mpqc_info MPQClcao_factory MPQCutil_keyval ${Boost_LIBRARIES} ${MADNESS_DISABLEPIE_LINKER_FLAG})

