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
  set(subdirectories
      ${subdirectories}
      madness)
  set(util_deplibs ${MADNESS_LIBRARIES})
endif()

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

set(util/OBJECTS ${objects} PARENT_SCOPE)
add_library(util ${objects})
install(TARGETS util DESTINATION lib)
target_link_libraries(util ${util_deplibs})
