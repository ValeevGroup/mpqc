set(sources
    decomposed_tensor_algebra.cpp
    decomposed_tensor_algebra.h
    decomposed_tensor.cpp
    decomposed_tensor.h
)

add_library(tensor ${sources})
target_link_libraries(tensor ${TA_LIBRARIES})
