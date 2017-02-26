set(sources
  group.cpp
  group.h
  petite_list.cpp
  petite_list.h
)

add_mpqc_library(math_groups sources sources "MPQCmath_eigen" "mpqc/math/groups")
