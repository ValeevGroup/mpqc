include(AddLibrary)

# elemental libs

set(sources
    eigensolver.cpp
)

add_mpqc_object_library(elemental ${sources})