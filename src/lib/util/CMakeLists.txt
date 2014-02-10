include(AddLibrary)

SUBDIRLIST(Subdirs CMAKE_CURRENT_SOURCE_DIR)

FOREACH(subdir ${Subdirs})
    ADD_SUBDIRECTORY(${subdir})
ENDFOREACH()

set(utilobjects
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
  list(APPEND utilobjects madness)
endif()

add_mpqc_library(util ${utilobjects})
install(TARGETS util DESTINATION lib)

if (HAVE_MADNESS)
  target_link_libraries(util ${MADNESS_LIBRARIES})
endif()

if (HAVE_MPI)
  target_link_libraries(util ${MPI_LIBRARIES})
endif()

# N.B. dlopen dependencies
target_link_libraries(util ${CMAKE_DL_LIBS})
 