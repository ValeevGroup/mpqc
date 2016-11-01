set(sources
    decomposed_tensor_algebra.cpp
    decomposed_tensor_algebra.h
    decomposed_tensor.cpp
    decomposed_tensor.h
)

add_library(MPQCmath_clr ${sources})
target_link_libraries(MPQCmath_clr ${TA_LIBRARIES})
