set(sources
  assert.h
  bug.cpp
  bug.h
  observer.h
  pool.h
  print.cpp
  print.h
  provider.h
  string.h
  task.h
  time.h
  timestamp.h
)

add_mpqc_library(util_misc sources sources "MADworld;tbb;MPQCutil_c++;MPQCutil_keyval" "mpqc/util/misc")
target_link_libraries(MPQCutil_misc PUBLIC MPQC_Eigen)
