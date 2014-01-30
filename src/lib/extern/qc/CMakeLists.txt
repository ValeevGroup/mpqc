# extern_qc library

add_library(extern_qc
  qc_fnc.cc
)
target_link_libraries(extern_qc chemistry)
install(TARGETS extern_qc DESTINATION lib)

if (MPQC_UNITTEST)
  set_property(SOURCE test_qcfnc0.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(test_qcfnc0 test_qcfnc0.cc)
  target_link_libraries(test_qcfnc0
    math
    util
    chemistry
  )
  add_test(test_qcfnc0 test_qcfnc0)
endif()

