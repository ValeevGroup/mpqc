include(AddLibrary)

# madness library

IF(MPQC_HAS_MADNESS)
    add_mpqc_object_library(madness
      init.cc
      world.cc
    )
ENDIF(MPQC_HAS_MADNESS)

# tests
if (MPQC_UNITTEST)
  add_executable(madtest madtest.cc $<TARGET_OBJECTS:class>)
  target_link_libraries(madtest util)
  target_link_libraries(madtest ${ELEMENTAL_LIBRARIES})
  target_link_libraries(madtest ${LAPACK_LIBRARIES})
  target_link_libraries(madtest ${Boost_LIBRARIES})
  add_test(madtest madtest)
endif()

