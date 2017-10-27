set(sources
  forcelink.h
  keyval.h
  keyval.cpp
)

add_mpqc_library(util_keyval sources sources "MPQCutil_c++;MPQCutil_meta;MPQCutil_core" "mpqc/util/keyval")

target_include_directories(MPQCutil_keyval PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(MPQCutil_keyval PUBLIC ${Boost_LIBRARIES})
