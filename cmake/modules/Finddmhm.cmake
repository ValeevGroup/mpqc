if (DMHM_INCLUDE_DIRS)

  set(DMHM_FOUND TRUE)

else ()

  find_path(DMHM_INCLUDE_DIR NAMES dmhm.hpp
      PATHS
      ${DMHM_PATH}/include
    )
  set(DMHM_INCLUDE_DIR ${DMHM_INCLUDE_DIR})

  if(DMHM_INCLUDE_DIR)
    set(DMHM_LIB_DIR ${DMHM_PATH}/lib)
    set(DMHM_DEFAULT_COMPONENT_LIST dmhm)

    foreach(_component ${DMHM_DEFAULT_COMPONENT_LIST})
      list(FIND DMHM_DEFAULT_COMPONENT_LIST ${DMHM_LIB_DIR}/${_component} _comp_found)

      #if(_comp_found EQUAL -1)
      #  message(FATAL_ERROR "Invalid DMHM component: ${_component}")
      #endif()

      if(_comp_found)
        set(DMHM_${_component}_FOUND TRUE)
      endif()
    endforeach()

    set(DMHM_LIBRARIES "")
    foreach(_component ${DMHM_DEFAULT_COMPONENT_LIST})
      message(STATUS " component = ${_component}")
      if(DMHM_${_component}_FOUND)
        list(APPEND DMHM_LIBRARIES ${_component})
        message(STATUS " Made it here ${_component}")
      endif()
      message(STATUS " dmhm lib ${DMHM_LIBRARIES}")
    endforeach()

    set(DMHM_FOUND TRUE)    
    
  endif(DMHM_INCLUDE_DIR)

  mark_as_advanced(DMHM_INCLUDE_DIRS DMHM_LIBRARIES)
  message(STATUS " INCLUDE = ${DMHM_INCLUDE_DIR}")
  SET(DMHM_INCLUDE_DIRS ${DMHM_INCLUDE_DIR} CACHE PATH "The dmhm include path.")
  SET(DMHM_LIBRARIES ${DMHM_LIBRARIES} CACHE PATH "The dmhm library path.")

endif()

