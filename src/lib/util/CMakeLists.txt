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
  set(util_deplibs ${MADNESS_LIBRARIES})
endif()

add_mpqc_library(util ${subdirectories})
install(TARGETS util DESTINATION lib)

target_link_libraries(util ${util_deplibs})
# N.B. dlopen dependencies
target_link_libraries(util ${CMAKE_DL_LIBS})
 