set(integral_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/screen_base.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/schwarz_screen.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/schwarz_screen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/cached_shell_info.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/cached_shell_info.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qqr_screening.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qqr_screening.cpp
)
add_library(integrals ${integral_sources})
#set_target_properties(integrals PROPERTIES LINKER_LANGUAGE C++)

