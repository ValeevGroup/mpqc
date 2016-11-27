set(sources
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
)

add_mpqc_library(util_misc sources sources "MADworld;util_c++" "mpqc/util/misc")
