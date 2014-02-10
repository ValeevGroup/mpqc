include(AddLibrary)

# ci library

IF(MPQC_CI)
    set(ccfiles
      ci.cc
      compute.cc
    )
    add_mpqc_object_library(ci ${ccfiles})
ENDIF()

# tests
if (MPQC_UNITTEST)
 # no tests yet
endif()

