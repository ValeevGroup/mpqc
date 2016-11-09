set(sources
  assert.h
  bug.h
  bug.cc
  exception.h
  exception.cc
  exenv.cc
  exenv.h
  formio.cc
  formio.h
  print.cpp
  print.h
  string.h
  time.h
  type_traits.h
)

add_library(MPQCutil_misc ${sources})
target_link_libraries(MPQCutil_misc MADworld)
