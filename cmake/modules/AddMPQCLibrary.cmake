macro(add_mpqc_library _name _rawlist_source_files _rawlist_public_header_files _dep_mpqc_comp _include_dir)

  set(_source_files )
  foreach(file ${${_rawlist_source_files}})
    if (file MATCHES ".*\\.cpp$")
      list(APPEND _source_files ${file})
    endif()
  endforeach()
  set(_public_header_files )
  foreach(file ${${_rawlist_public_header_files}})
    if (NOT file MATCHES ".*\\.cpp$")
      list(APPEND _public_header_files ${file})
    endif()
  endforeach()
  
  #message (STATUS "in add_mpqc_library: name=${_name} source_files=${_source_files} public_header_files=${_public_header_files} deps=${_dep_mpqc_comp} include_dir=${_include_dir}")
  if (NOT _source_files) # no sources given?
    message (FATAL_ERROR "add_mpqc_library: no sources given; probably want add_mpqc_hdr_library instead")
  endif()

  set(_libname MPQC${_name})
  
  add_library(${_libname}-obj OBJECT ${_source_files})
  #cotire(${_libname}-obj)
  add_library(${_libname} $<TARGET_OBJECTS:${_libname}-obj>)
  if(BUILD_SHARED_LIBS)
    set_target_properties(${_libname}-obj PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  endif()

  target_include_directories(${_libname} PUBLIC
    $<INSTALL_INTERFACE:${MPQC_INSTALL_INCLUDEDIR}>)
  
  # Pass the private ${_libname} compile flags to ${_libname}-obj  
  target_compile_definitions(${_libname}-obj PRIVATE 
      $<TARGET_PROPERTY:${_libname},COMPILE_DEFINITIONS>)
  target_include_directories(${_libname}-obj PRIVATE 
      $<TARGET_PROPERTY:${_libname},INCLUDE_DIRECTORIES>)
  target_compile_options(${_libname}-obj PRIVATE 
      $<TARGET_PROPERTY:${_libname},COMPILE_OPTIONS>)
  
  # this does not work with hierarchies of header files
  # see here for possible workaround if frameworks are really needed:
  #     http://cmake.3232098.n2.nabble.com/Install-header-directory-hierarchy-td5638507.html
  # set_target_properties(${_libname} PROPERTIES PUBLIC_HEADER "${_public_header_files}")
  # install manually
  foreach ( file ${_public_header_files} )
    get_filename_component( dir ${file} DIRECTORY )
    install( FILES ${file} DESTINATION ${MPQC_INSTALL_INCLUDEDIR}/${_include_dir}/${dir} COMPONENT ${_name})
  endforeach()
  
  
  # Add library to the list of installed components
  install(TARGETS ${_libname} EXPORT mpqc
      COMPONENT ${_name}
      LIBRARY DESTINATION "${MPQC_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${MPQC_INSTALL_LIBDIR}"
      INCLUDES DESTINATION "${MPQC_INSTALL_INCLUDEDIR}")
  
  # Create a target to install the component
  add_custom_target(install-${_name}
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${_name} -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing ${_name} library components")
  add_dependencies(install-${_name} ${_libname})

  set(LINK_FLAGS "")
  foreach(_dep ${_dep_mpqc_comp})
    # for MPQC dependencies make sure they have been processed already!
    if (${_dep} MATCHES "^MPQC" AND NOT DEFINED ${_dep}_is_mpqc_hdr_lib)
      message(FATAL_ERROR "add_mpqc_library given ${_dep} as dependence of ${_name}, but not defined yet!")
    endif()
    
    if(TARGET install-${_dep})
      add_dependencies(install-${_name} install-${_dep})
    endif()
    if(TARGET ${_dep})
      target_compile_definitions(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_DEFINITIONS>)
      target_include_directories(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_INCLUDE_DIRECTORIES>)
      target_compile_options(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_OPTIONS>)
      if (${_dep}_is_mpqc_hdr_lib)
        target_link_libraries(${_libname} INTERFACE ${_dep})
      else()
        target_link_libraries(${_libname} PUBLIC ${_dep})
      endif()      

      # import LINK_FLAGS from dependent, unless it's an interface library
      get_property(${_dep}_target_type TARGET ${_dep} PROPERTY TYPE SET)
      if (NOT ${_dep}_is_mpqc_hdr_lib AND NOT ${_dep}_target_type STREQUAL "INTERFACE_LIBRARY")
        get_property(_dep_LINK_FLAGS_SET TARGET ${_dep} PROPERTY LINK_FLAGS SET)
        if (_dep_LINK_FLAGS_SET)
          get_property(_dep_LINK_FLAGS TARGET ${_dep} PROPERTY LINK_FLAGS)
          set(LINK_FLAGS "${LINK_FLAGS} ${_dep_LINK_FLAGS}")
        endif ()
      endif ()
      
    endif()
  endforeach()
  set_target_properties(${_libname} PROPERTIES LINK_FLAGS "${LINK_FLAGS}") 
 
  # Use current CMAKE_CXX_FLAGS to compile targets dependent on this library
  string (REPLACE " " ";" CMAKE_CXX_FLAG_LIST "${CMAKE_CXX_FLAGS}")
  target_compile_options(${_libname} INTERFACE $<INSTALL_INTERFACE:${CMAKE_CXX_FLAG_LIST}>)

  set(${_libname}_is_mpqc_hdr_lib FALSE CACHE INTERNAL "")
endmacro()

macro(add_mpqc_hdr_library _name _public_header_files _dep_mpqc_comp _include_dir)

  set(_libname MPQC${_name})

  # make INTERFACE library
  add_library(${_libname} INTERFACE)
  
  target_include_directories(${_libname} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${MPQC_INSTALL_INCLUDEDIR}>
  )

  # PUBLIC_HEADER property is not supported on INTERFACE targets, and do not work
  # with hierarchies of headers
  # hence have to do install manually
  foreach ( file ${${_public_header_files}} )
    get_filename_component( dir ${file} DIRECTORY )
    install( FILES ${file} DESTINATION ${MPQC_INSTALL_INCLUDEDIR}/${_include_dir}/${dir} COMPONENT ${_name})
  endforeach()
  
  # Add library to the list of installed components
  install(TARGETS ${_libname} EXPORT mpqc
          COMPONENT ${_name})
  
  # Create a target to install the component
  add_custom_target(install-${_name}
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${_name} -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing ${_name} library components")
  add_dependencies(install-${_name} ${_libname})

  foreach(_dep ${_dep_mpqc_comp})
    if(TARGET install-${_dep})
      add_dependencies(install-${_name} install-${_dep})
    endif()
    if(TARGET ${_dep})
        target_compile_definitions(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_DEFINITIONS>)
        target_include_directories(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_INCLUDE_DIRECTORIES>)
        target_compile_options(${_libname} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_OPTIONS>)
      if (${_dep}_is_mpqc_hdr_lib)
        target_link_libraries(${_libname} INTERFACE ${_dep})
      else()
        target_link_libraries(${_libname} PUBLIC ${_dep})
      endif()
    endif()
  endforeach()
  
  # Use current CMAKE_CXX_FLAGS to compile targets dependent on this library
  string (REPLACE " " ";" CMAKE_CXX_FLAG_LIST "${CMAKE_CXX_FLAGS}")
  target_compile_options(${_libname} INTERFACE $<INSTALL_INTERFACE:${CMAKE_CXX_FLAG_LIST}>)
  
  set(${_libname}_is_mpqc_hdr_lib TRUE CACHE INTERNAL "")
endmacro()
