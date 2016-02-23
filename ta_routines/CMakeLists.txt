set(sources 
    cholesky_inverse.h
    cholesky_inverse.cpp
        )

add_library(ta_routines ${sources})
target_link_libraries(ta_routines ${TA_LIBRARIES})
