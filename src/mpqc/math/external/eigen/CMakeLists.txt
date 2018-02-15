set(sources
  eigen.h
  util.h
  util.cpp
)

add_mpqc_hdr_library(math_eigen sources "" "mpqc/math/external/eigen")
