set(sources
  keyval.hpp
  keyval.cpp
)

add_library(MPQCkeyval ${sources})
target_link_libraries(MPQCkeyval ${BOOST_LIBRARIES})

#  set_property(SOURCE keyvaltest.cpp PROPERTY COMPILE_DEFINITIONS
#      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#  add_executable(keyvaltest keyvaltest.cpp)
#  target_link_libraries(keyvaltest
#    MPQCkeyval
#  )
#  add_test(keyvaltest keyvaltest)
