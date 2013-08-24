include(AddLibrary)

# distarray4 library

set(sources
  distarray4.cc
  distarray4_memgrp.cc
  distarray4_node0file.cc
)

if (HAVE_MPI)
  list(APPEND sources distarray4_mpiiofile.cc)
endif()

add_mpqc_object_library(distarray4 ${sources})

# tests

