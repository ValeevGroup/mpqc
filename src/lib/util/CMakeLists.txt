include(AddLibrary)

set(subdirectories
  class
  container
  group
  keyval
  misc
  options
  ref
  render
  state
)

if (HAVE_MADNESS)
  list(APPEND subdirectories madness)
endif()

add_mpqc_library(util ${subdirectories})
install(TARGETS util DESTINATION lib)

if (HAVE_MADNESS)
  target_link_libraries(util ${MADNESS_LIBRARIES})
endif()

if (HAVE_MPI)
  target_link_libraries(util ${MPI_LIBRARIES})
endif()

# N.B. dlopen dependencies
target_link_libraries(util ${CMAKE_DL_LIBS})
 