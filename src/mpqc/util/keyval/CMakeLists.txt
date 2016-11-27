set(sources
  forcelink.h
  keyval.h
  keyval.cpp
)

add_mpqc_library(util_keyval sources sources "util_c++" "mpqc/util/keyval")

target_include_directories(MPQCutil_keyval PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(MPQCutil_keyval PUBLIC ${Boost_LIBRARIES})
