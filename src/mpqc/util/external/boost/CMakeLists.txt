set(sources
  dummy.cpp
  small_vector.h
)

add_mpqc_library(util_boost sources sources "MADworld" "mpqc/util/external/boost")
target_include_directories(MPQCutil_boost PUBLIC ${Boost_INCLUDE_DIRS})
