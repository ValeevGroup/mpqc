set(sources
  basis.h
  libintshell.h
  clusterable_shell.h
)

add_executable(libintshell_test libintshell_test.cpp)
target_link_libraries(libintshell_test ${LIBINT_LIBRARIES})
