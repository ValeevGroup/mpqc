set(sources
        linkage.cpp
        phf.h
        phf.cpp
        )

add_library(MPQCphf ${sources})
target_link_libraries(MPQCphf MPQCwfn MPQCmath_ta MPQCmath_linalg MPQCmath_clr)
