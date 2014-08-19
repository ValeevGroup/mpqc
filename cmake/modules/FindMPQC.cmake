if (MPQC_INCLUDE_DIRS)

  set(MPQC_FOUND TRUE)

else (MPQC_INCLUDE_DIRS)

  set (MPQC_FOUND "NO")

  if(MPQC_INSTALL_DIR)
    set(_MPQC_INSTALL_DIR ${MPQC_INSTALL_DIR})
  endif(MPQC_INSTALL_DIR)

  set(MPQC_INC_SEARCH_DIR ${_MPQC_INSTALL_DIR}/include)

  find_path(MPQC_INCLUDE_DIR NAME mpqc/mpqc_config.h
      HINTS
      ${MPQC_INC_SEARCH_DIR}
    )

  mark_as_advanced(MPQC_INCLUDE_DIR)

  set(_MPQC_LIB_NAME "chemistry")
  set(_MPQC_LIBRARY_DIR ${_MPQC_INSTALL_DIR}/lib)

  find_library(MPQC_LIBRARY ${_MPQC_LIB_NAME} HINTS ${_MPQC_LIBRARY_DIR})
  get_filename_component(MPQC_LIBRARY_DIR ${MPQC_LIBRARY} PATH)

  mark_as_advanced(MPQC_LIBRARY)

  if(MPQC_INCLUDE_DIR)
    if(MPQC_LIBRARY)
      set (MPQC_FOUND "YES")
      set (MPQC_LIBRARIES ${MPQC_LIBRARY} CACHE PATH "MPQC libraries")
      set (MPQC_INCLUDE_DIRS ${MPQC_INCLUDE_DIR} CACHE PATH "MPQC include directory" FORCE)
      set (MPQC_LIBRARY_DIRS ${MPQC_LIBRARY_DIR} CACHE PATH "MPQC library directory" FORCE)
      mark_as_advanced(MPQC_INCLUDE_DIRS MPQC_LIBRARY_DIRS MPQC_LIBRARIES)
      message(STATUS "Found MPQC libraries ${MPQC_LIBRARIES}")
    endif(MPQC_LIBRARY)
  endif(MPQC_INCLUDE_DIR)

  if(NOT MPQC_FOUND)
    message("ERROR: MPQC NOT found!")
    message(STATUS "Looked for MPQC in ${_MPQC_INSTALL_DIR}")
      if(MPQC_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find MPQC library.")
      endif(MPQC_FIND_REQUIRED)
  endif(NOT MPQC_FOUND)

endif(MPQC_INCLUDE_DIRS)


