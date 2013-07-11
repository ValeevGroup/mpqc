# sources

add_subdirectory(mpqc)

if (HAVE_LIBINT2)
  add_subdirectory(pt2r12)
endif()
