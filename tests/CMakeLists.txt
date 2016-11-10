ADD_SUBDIRECTORY(unit)
ADD_SUBDIRECTORY(validation)

add_custom_target(check
                  COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target check_validation
                  COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target check_unit
                 )
