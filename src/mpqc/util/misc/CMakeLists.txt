set(sources
  exception.hpp
  exception.cpp
)

add_library(MPQCmisc ${sources})
target_link_libraries(MPQCmisc ${BOOST_LIBRARIES})
