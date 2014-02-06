include(AddLibrary)

# ci library

IF(MPQC_CI)
    add_mpqc_object_library(ci
      ci.cc
      compute.cc
    )
ENDIF()

# tests
if (MPQC_UNITTEST)
 # no tests yet
endif()

