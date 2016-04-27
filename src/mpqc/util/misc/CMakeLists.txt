set(sources
  exception.hpp
  type_traits.h
  exception.cpp
)

add_library(MPQCmisc ${sources})
target_link_libraries(MPQCmisc ${BOOST_LIBRARIES})
