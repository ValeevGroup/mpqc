set(sources
    energy.h
    energy.cpp
)

add_library(MPQCproperties ${sources})
target_link_libraries(MPQCproperties MPQCwfn)
