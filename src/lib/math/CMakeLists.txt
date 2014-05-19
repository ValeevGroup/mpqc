include(AddLibrary)
  
SUBDIRLIST(Subdirs CMAKE_CURRENT_SOURCE_DIR)

FOREACH(subdir ${Subdirs})
    ADD_SUBDIRECTORY(${subdir})
ENDFOREACH()

set(mathobjects
  distarray4
  isosurf
  mmisc
  optimize
  scmat
  symmetry
  )
  
if(MPQC_NEW_FEATURES)  
  LIST(APPEND mathobjects elemental)
endif(MPQC_NEW_FEATURES)
  
add_mpqc_library(math ${mathobjects})

add_dependencies(math util)
target_link_libraries(math util ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

if(MPQC_NEW_FEATURES)
    TARGET_LINK_LIBRARIES(math ${ELEMENTAL_LIBRARIES})
endif(MPQC_NEW_FEATURES)

set(math/OBJECTS ${objects} PARENT_SCOPE)
install(TARGETS math DESTINATION lib)
