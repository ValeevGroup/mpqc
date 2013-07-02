
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h
)


# set(CC_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/cc/integrals.cpp")
# #file(GLOB CC_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/cc/*.cpp")
# add_library(mpqc/cc OBJECT
#   ${CC_SOURCES}
# )

# add_executable(test-mpqc EXCLUDE_FROM_ALL test.cpp)

# target_link_libraries(test-mpqc ${Boost_LIBRARIES})
# #target_link_libraries(test-mpqc ${LAPACK_LIBRARIES})
# target_link_libraries(test-mpqc ${PYTHON_LIBRARIES})
# target_link_libraries(test-mpqc ${LIBRARIES})
# target_link_libraries(test-mpqc ${CMAKE_DL_LIBS})
