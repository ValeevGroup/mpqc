include(AddLibrary)

# ci library

add_mpqc_object_library(ci
  ci.cc
  compute.cc
)

# tests
if (MPQC_UNITTEST)
 # no tests yet
endif()

