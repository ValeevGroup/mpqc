set(sources
  array_info.cpp
  array_info.h
  tensor_store.h
  util.h
  util.cpp
)

add_mpqc_library(math_ta sources sources "tiledarray" "mpqc/math/external/tiledarray")
