set(sources
    array_info.h
    array_info.cpp
    memory_usage.h
    memory_usage.cpp
    json_handling.h
    json_handling.cpp
)

add_library(utility ${sources})
target_link_libraries(utility ${TA_LIBRARIES})

