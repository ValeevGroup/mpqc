set(sources
  basisset.h
  atom_basisset.h
  basisset.cpp
  clusterable_shell.h
)

add_executable(libintshell_test libintshell_test.cpp)
target_link_libraries(libintshell_test ${LIBINT_LIBRARIES})
