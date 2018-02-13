set(sources
  eigen.h
  util.h
  util.cpp
)

add_mpqc_library(math_eigen sources sources "" "mpqc/math/external/eigen")
target_link_libraries(MPQCmath_eigen PUBLIC MPQC_Eigen)