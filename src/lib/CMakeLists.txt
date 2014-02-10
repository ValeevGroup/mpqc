add_subdirectory(mpqc)
add_subdirectory(chemistry)
add_subdirectory(math)
add_subdirectory(util)
add_subdirectory(extern)

# library headers
file(GLOB_RECURSE headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} chemistry/*.h math/*.h util/*h)
foreach (_file ${headers})
  get_filename_component(_path ${_file} PATH)
  #message(STATUS "gonna installl ${_file} in include/mpqc/${_path}" )
  install(FILES ${_file} DESTINATION include/mpqc/${_path})
endforeach()
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mpqc_config.h DESTINATION include/mpqc)