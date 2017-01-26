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
  observer.h
  pool.h
  print.cpp
  print.h
  string.h
  time.h
  timestamp.h
)

add_mpqc_library(util_misc sources sources "MADworld;util_c++" "mpqc/util/misc")
target_link_libraries(MPQCutil_misc PUBLIC ${TBB_LIBRARIES})
