set(sources
  print.cpp
  print.h
  string.h
  time.h
  type_traits.h
)

add_library(MPQCutil_misc ${sources})
target_link_libraries(MPQCutil_misc MADworld)
