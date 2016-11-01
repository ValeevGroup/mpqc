set(sources
  json_handling.cpp
  json_handling.h
  string.h
  time.h
  type_traits.h
)

add_library(MPQCmisc ${sources})
target_link_libraries(MPQCmisc MADworld)
