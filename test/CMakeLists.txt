cmake_minimum_required(VERSION 2.8)

# skeleton for CTest - implement in future
# #add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
# enable_testing()

# file(GLOB REF RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ref/*.in)

# foreach(input ${REF})
#   add_test(${input} ../mpqc.py "-X ${CMAKE_CURRENT_SOURCE_DIR}/../src/bin/mpqc/mpqc" ${input})
# endforeach()

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/check.mk.in
  ${CMAKE_CURRENT_BINARY_DIR}/check.mk
  @ONLY
)

# make.sh runs tests for every integral factory
set(MPQC_INTEGRAL_FACTORY_LIST "intv3")
if (HAVE_LIBINT1)
  list(APPEND MPQC_INTEGRAL_FACTORY_LIST "cints")
endif (HAVE_LIBINT1)
if (HAVE_LIBINT2)
  list(APPEND MPQC_INTEGRAL_FACTORY_LIST "libint2")
endif (HAVE_LIBINT2)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/check.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/check.cmake
  @ONLY
)

add_custom_target(
  check0
  COMMAND ${CMAKE_COMMAND} -DCHECKLEVEL=check0 -P ${CMAKE_CURRENT_BINARY_DIR}/check.cmake
)

add_custom_target(
  check1
  COMMAND ${CMAKE_COMMAND} -DCHECKLEVEL=check1 -P ${CMAKE_CURRENT_BINARY_DIR}/check.cmake
)

add_custom_target(
  check2
  COMMAND ${CMAKE_COMMAND} -DCHECKLEVEL=check2 -P ${CMAKE_CURRENT_BINARY_DIR}/check.cmake
)

add_custom_target(
  check
  DEPENDS check0
)

# checking
add_custom_target(
  check_clean
  COMMAND make -f check.mk check_clean
)
