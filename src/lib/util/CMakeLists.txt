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

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

set(util/OBJECTS ${objects} PARENT_SCOPE)
add_library(util ${objects})
