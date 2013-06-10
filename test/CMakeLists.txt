cmake_minimum_required(VERSION 2.8)

# skeleton for CTest - implement in future
# #add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
# enable_testing()

# file(GLOB REF RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ref/*.in)

# foreach(input ${REF})
#   add_test(${input} ../mpqc.py "-X ${CMAKE_CURRENT_SOURCE_DIR}/../src/bin/mpqc/mpqc" ${input})
# endforeach()


add_custom_target(
  check0
  COMMAND make -f check.mk check0
)

add_custom_target(
  check1
  COMMAND make -f check.mk check1
)

add_custom_target(
  check2
  COMMAND make -f check.mk check2
)

add_custom_target(
  check
  DEPENDS check0 check1 check2
)

# checking
add_custom_target(
  check_clean
  COMMAND make -f check.mk check_clean
)
