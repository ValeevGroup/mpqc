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

if (HAVE_ELEMENTAL)
  list(APPEND utilobjects elemental)
endif()

add_mpqc_library(util ${utilobjects})
install(TARGETS util DESTINATION lib)

target_link_libraries(util ${MADNESS_LIBRARIES})
target_link_libraries(util ${ELEMENTAL_LIBRARIES})
target_link_libraries(util ${MPI_LIBRARIES})
target_link_libraries(util ${Boost_LIBRARIES})
target_link_libraries(util ${ARMCI_LIBRARIES})
target_link_libraries(util ${Threads_LIBRARIES})

# N.B. dlopen dependencies
target_link_libraries(util ${CMAKE_DL_LIBS})
 