include(AddLibrary)

set(subdirectories
  distarray4
  isosurf
  mmisc
  optimize
  scmat
  symmetry
  )

add_mpqc_library(math ${subdirectories})

add_dependencies(math util)
target_link_libraries(math util ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
set(math/OBJECTS ${objects} PARENT_SCOPE)
install(TARGETS math DESTINATION lib)
