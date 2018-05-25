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

add_mpqc_library(util_misc sources sources "MADworld;MPQCutil_c++;MPQCutil_keyval" "mpqc/util/misc")
# since FindTBB does not create an interface target, we can't just use it as a dependence in add_mpqc_library call
target_link_libraries(MPQCutil_misc PUBLIC tbb MPQC_Eigen)
