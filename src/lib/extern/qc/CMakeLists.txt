# extern_qc library

add_library(extern_qc
  qc_fnc.cc
)
target_link_libraries(extern_qc chemistry)
install(TARGETS extern_qc DESTINATION lib)

set_property(SOURCE test_qcfnc0.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(test_qcfnc0 EXCLUDE_FROM_ALL test_qcfnc0.cc)
add_test(test_qcfnc0 test_qcfnc0)

