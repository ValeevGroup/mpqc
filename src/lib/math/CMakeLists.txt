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
add_dependencies(math util)
target_link_libraries(math util ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
set(math/OBJECTS ${objects} PARENT_SCOPE)
install(TARGETS math DESTINATION lib)
