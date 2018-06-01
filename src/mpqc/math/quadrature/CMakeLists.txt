set(sources 
 	gaussian.h
 	gaussian.cpp
)

add_mpqc_library(math_quadrature sources sources "MPQCmath_eigen" "mpqc/math/quadrature")
