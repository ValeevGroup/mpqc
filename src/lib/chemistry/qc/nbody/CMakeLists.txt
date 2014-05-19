include(AddLibrary)

# nbody library

SET(ccfiles
  nbwfn.cc
  ref.cc
)

add_mpqc_object_library(nbody ${ccfiles})

# tests

