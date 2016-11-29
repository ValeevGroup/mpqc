set(sources
  parallel_break_point.h
  parallel_file.h
  parallel_file.cpp
  parallel_print.h
)

add_mpqc_library(util_mad sources sources "MADworld" "mpqc/util/external/madworld")
