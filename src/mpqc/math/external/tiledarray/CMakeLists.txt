set(sources
  array_info.h
  array_info.cpp
  tensor_store.h
)

add_library(MPQCmath_ta ${sources})
target_link_libraries(MPQCmath_ta tiledarray)
