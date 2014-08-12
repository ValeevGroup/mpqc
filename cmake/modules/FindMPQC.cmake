if (MPQC_INCLUDE_DIRS)

  set(MPQC_FOUND TRUE)

else ()

  SET(MPQC_INCLUDE_DIR "")
  find_path(MPQC_INCLUDE_DIR NAMES mpqc_config.h
      PATHS
      ${MPQC_PATH}
      ${MPQC_PATH}/include
      ${MPQC_PATH}/include/mpqc
      PATH_SUFFIXES mpqc
    )
  message(STATUS "MPQC_INCLUDE_DIR = ${MPQC_INCLUDE_DIR}")

  if(MPQC_INCLUDE_DIR)
    set(MPQC_LIB_DIR ${MPQC_PATH}/lib)
    set(MPQC_DEFAULT_COMPONENT_LIST chemistry)

    foreach(_component ${MPQC_DEFAULT_COMPONENT_LIST})
      list(FIND MPQC_DEFAULT_COMPONENT_LIST ${MPQC_LIB_DIR}/${_component} _comp_found)

      #if(_comp_found EQUAL -1)
      #  message(FATAL_ERROR "Invalid DMHM component: ${_component}")
      #endif()

      if(_comp_found)
        set(MPQC_${_component}_FOUND TRUE)
      endif()
    endforeach()

    foreach(_component ${MPQC_DEFAULT_COMPONENT_LIST})
      if(MPQC_${_component}_FOUND)
        list(APPEND MPQC_LIBRARIES ${_component})
      endif()
    endforeach()

    set(MPQC_FOUND TRUE)    
    
  endif(MPQC_INCLUDE_DIR)

  mark_as_advanced(MPQC_INCLUDE_DIR MPQC_LIBRARIES)
  SET(MPQC_INCLUDE_DIRS ${MPQC_INCLUDE_DIR} CACHE PATH "The mpqc include path.")
  SET(MPQC_LIBRARIES ${MPQC_LIBRARIES} CACHE PATH "The mpqc library path.")

endif()

