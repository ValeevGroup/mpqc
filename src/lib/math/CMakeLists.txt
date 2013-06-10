set(subdirectories
  distarray4
  isosurf
  mmisc
  optimize
  scmat
  symmetry
  )

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

add_library(math ${objects})
set(math/OBJECTS ${objects} PARENT_SCOPE)
