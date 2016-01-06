set(sources
    utility.cpp
    utility.h
)

add_library(f12 ${sources})
target_link_libraries(f12 ${TA_LIBRARIES})
