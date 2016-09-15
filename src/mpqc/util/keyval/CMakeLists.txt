set(sources
  keyval.hpp
  keyval.cpp
)

add_library(MPQCkeyval ${sources})
target_link_libraries(MPQCkeyval ${BOOST_LIBRARIES})
