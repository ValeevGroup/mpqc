include(AddLibrary)

# misc library

set(ccfiles
  bug.cc
  compute.cc
  consumableresources.cc
  exception.cc
  exenv.cc
  formio.cc
  ieee.cc
  regtime.cc
  runnable.cc
  scexception.cc
  units.cc
)
if (MPQC_NEW_FEATURES)
  list(APPEND ccfiles
    thread_timer.cc
    xml.cc
    xmlwriter.cc
  )
endif (MPQC_NEW_FEATURES)

add_mpqc_object_library(misc ${ccfiles})

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE autovectest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(autovectest autovectest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(autovectest
    util
    ${LINK_LIBRARIES}
  )
  add_test(autovectest autovectest)

  set_property(SOURCE bugtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(bugtest bugtest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(bugtest
    util
  )
  add_test(bugtest bugtest)

  set_property(SOURCE comptest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(comptest comptest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(comptest
    util
  )
  add_test(comptest comptest)

  set_property(SOURCE formiotest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(formiotest formiotest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(formiotest
    util
  )
  add_test(formiotest formiotest)

  set_property(SOURCE regtimetest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(regtimetest regtimetest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(regtimetest
    util
  )
  add_test(regtimetest regtimetest)

  set_property(SOURCE scinttest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(scinttest scinttest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(scinttest
    util
  )
  add_test(scinttest scinttest)

  set_property(SOURCE unittest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(unittest unittest.cc $<TARGET_OBJECTS:misc>)
  target_link_libraries(unittest
    util
  )
  add_test(unittest unittest)

endif()
