set(sources
  assert.h
  bug.cpp
  bug.h
  exception.cpp
  exception.h
  exenv.cpp
  exenv.h
  formio.cpp
  formio.h
  print.cpp
  print.h
  string.h
  time.h
  type_traits.h
)

add_library(MPQCutil_misc ${sources})
target_link_libraries(MPQCutil_misc MADworld)
