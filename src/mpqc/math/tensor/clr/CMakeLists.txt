set(sources
array_to_eigen.h
decomposed_tensor_addition.h
decomposed_tensor_algebra.cpp
decomposed_tensor_algebra.h
decomposed_tensor_gemm_helper.h
decomposed_tensor_gemm.h
decomposed_tensor_multiplication.h
decomposed_tensor_nonintrusive_interface.h
decomposed_tensor_subtraction.h
decomposed_tensor_unary.h
decomposed_tensor.cpp
decomposed_tensor.h
minimize_storage.h
tensor_fwd.h
tensor_train.h
tensor_transforms.h
tile_algebra.h
tile.h
vector_localization.h
cp_als.h)

add_mpqc_library(math_clr sources sources "MPQCmath_eigen;MPQCmath_lapack" "mpqc/math/tensor/clr")
