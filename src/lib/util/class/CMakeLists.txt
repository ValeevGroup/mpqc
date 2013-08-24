include(AddLibrary)

# class library
add_mpqc_object_library(class class.cc)

# tests

set_property(SOURCE classtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(classtest EXCLUDE_FROM_ALL classtest.cc $<TARGET_OBJECTS:class>)
add_test(classtest classtest)

# set_property(SOURCE scextest.cc PROPERTY COMPILE_DEFINITIONS
#     -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
# add_executable(scextest EXCLUDE_FROM_ALL scextest.cc $<TARGET_OBJECTS:class>)
# add_test(scextest scextest)

