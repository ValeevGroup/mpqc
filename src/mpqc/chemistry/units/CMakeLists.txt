set(sources
        units.h
        units.cpp
)

add_library(MPQCunits ${sources})
target_link_libraries(MPQCunits MPQCutil_misc)
