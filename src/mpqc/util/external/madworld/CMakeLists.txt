set(sources
  parallel_break_point.h
  parallel_file.h
  parallel_file.cpp
  parallel_print.h
)

add_library(MPQCutil_mad ${sources})
target_link_libraries(MPQCutil_mad ${TA_LIBRARIES})
