set(sources
    decomposed_tensor_algebra.cpp
    decomposed_tensor_algebra.h
    decomposed_tensor.cpp
    decomposed_tensor.h
)

add_library(MPQCmath_tensor_clr ${sources})
target_link_libraries(MPQCmath_tensor_clr ${TA_LIBRARIES})
