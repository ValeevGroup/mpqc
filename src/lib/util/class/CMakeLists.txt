include(AddLibrary)

# class library
add_mpqc_object_library(class class.cc)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE classtest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(classtest classtest.cc $<TARGET_OBJECTS:class>)
  target_link_libraries(classtest util)
  add_test(classtest classtest)
endif()

