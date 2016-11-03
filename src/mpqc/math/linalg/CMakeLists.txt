set(sources 
    cholesky_inverse.h
    cholesky_inverse.cpp
)

add_library(MPQCmath_linalg ${sources})
target_link_libraries(MPQCmath_linalg ${TA_LIBRARIES})
