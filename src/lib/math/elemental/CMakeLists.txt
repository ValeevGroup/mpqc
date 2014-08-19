include(AddLibrary)

# elemental libs

IF(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)
  set(sources
      eigensolver.cpp
  )
ENDIF(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)

add_mpqc_object_library(elemental ${sources})
