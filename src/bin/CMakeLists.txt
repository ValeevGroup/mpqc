# sources

add_subdirectory(mpqc)
add_subdirectory(TA_Ints_test)

if (MPQC_R12)
  add_subdirectory(pt2r12)
endif()
