set(sources 
    cholesky_inverse.h
    conditioned_orthogonalizer.h
    diagonal_array.h
    eigen_value_estimation.h
    inverse.h
    sqrt_inv.h
)

add_mpqc_hdr_library(math_linalg sources sources "mpqc/math/linalg")
