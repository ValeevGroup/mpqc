set(integral_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/screen_base.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/schwarz_screen.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/schwarz_screen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/cached_shell_info.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/cached_shell_info.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qqr_screening.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qqr_screening.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qvl_screening.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qvl_screening.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qvl_shell_info.h
    ${CMAKE_CURRENT_SOURCE_DIR}/screening/qvl_shell_info.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/task_integral_kernels.cpp
)

add_library(integrals ${integral_sources})
target_link_libraries(integrals ${TA_LIBRARIES})

