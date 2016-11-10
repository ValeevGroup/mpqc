set(sources
  keyval.h
  keyval.cpp
)

add_library(MPQCutil_keyval ${sources})
target_link_libraries(MPQCutil_keyval ${Boost_LIBRARIES})
