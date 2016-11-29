set(sources 
    cholesky_inverse.h
    cholesky_inverse.cpp
    conditioned_orthogonalizer.h
    diagonal_array.h
    eigen_value_estimation.h
    inverse.h
    sqrt_inv.h
)

add_mpqc_library(math_linalg sources sources "MPQCmath_eigen;MPQCmath_ta" "mpqc/math/linalg")
