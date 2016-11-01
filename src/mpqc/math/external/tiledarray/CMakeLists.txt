set(sources
  array_info.cpp
  array_info.h
  tensor_store.h
)

add_library(MPQCmath_ta ${sources})
target_link_libraries(MPQCmath_ta ${TA_LIBRARIES})
