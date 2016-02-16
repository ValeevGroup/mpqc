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
string(REPLACE ";" " " MPQC_INTEGRAL_FACTORY_STRLIST "${MPQC_INTEGRAL_FACTORY_LIST}")
message (STATUS "MPQC_INTEGRAL_FACTORY_STRLIST=${MPQC_INTEGRAL_FACTORY_STRLIST}")
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/make.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/make.sh
  @ONLY
)

add_custom_target(
  check0
  COMMAND ./make.sh check0
)

add_custom_target(
  check1
  COMMAND ./make.sh check1
)

add_custom_target(
  check2
  COMMAND ./make.sh check2
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
