# extern_qc library

add_library(extern_qc
  qc_fnc.cc
)
target_link_libraries(extern_qc chemistry)
install(TARGETS extern_qc DESTINATION lib)

if (MPQC_UNITTEST)
  add_executable(test_qcfnc0 test_qcfnc0.cc qc_fnc.cc)
  target_link_libraries(test_qcfnc0
    chemistry
  )
  add_test(test_qcfnc0 test_qcfnc0)
endif()

