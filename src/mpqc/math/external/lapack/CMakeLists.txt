set(sources
  lapack.h
  lapack.cpp
)

add_mpqc_library(math_lapack sources sources "tiledarray" "mpqc/math/external/lapack")
