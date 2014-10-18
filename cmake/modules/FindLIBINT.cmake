if (LIBINT_INCLUDE_DIRS)
  set (LIBINT_FOUND TRUE)

else (LIBINT_INCLUDE_DIRS)
  set (LIBINT_FOUND FALSE)

  if(LIBINT_INSTALL_DIR)
    set(_LIBINT_INSTALL_DIR ${LIBINT_INSTALL_DIR})
  endif(LIBINT_INSTALL_DIR)

  set(LIBINT_INC_SEARCH_DIR ${_LIBINT_INSTALL_DIR}/include)
  message(STATUS "Searching for includes in ${LIBINT_INC_SEARCH_DIR}")

  find_path(LIBINT_INCLUDE_DIR NAMES libint2.h
    HINTS
    ${LIBINT_INC_SEARCH_DIR})

  find_path(LIBINT_INCLUDE_LIBINT2_DIR NAMES cxxapi.h
    HINTS
    ${LIBINT_INC_SEARCH_DIR}/libint2)

  set(LIBINT_INCLUDE_DIR ${LIBINT_INCLUDE_DIR} ${LIBINT_INCLUDE_LIBINT2_DIR})

  message(STATUS "LIBINT INCLUDE DIRS = ${LIBINT_INCLUDE_DIR}")

  mark_as_advanced(LIBINT_INCLUDE_DIR)

  set(_LIBINT_LIB_NAMES
    "int2"
    )

  set(_LIBINT_LIBRARY_DIR ${_LIBINT_INSTALL_DIR}/lib)
  set(LIBINT_LIBRARY "")

  foreach(_lib ${_LIBINT_LIB_NAMES})
    set(current_lib NOTFOUND)
    find_library(current_lib ${_lib} HINTS ${_LIBINT_LIBRARY_DIR})
    if(NOT current_lib)
      message(FATAL_ERROR "TileClusterChem could not find Libint lib: ${_lib}")
    endif(NOT current_lib)
    message(STATUS " Found Libint lib: ${current_lib}")
    LIST(APPEND LIBINT_LIBRARY ${current_lib})
  endforeach()

  set(LIBINT_LIBRARY_DIR ${_LIBINT_LIBRARY_DIR})
  mark_as_advanced(LIBINT_LIBRARY)

  if(LIBINT_INCLUDE_DIR)
    if(LIBINT_LIBRARY)
      set(LIBINT_FOUND TRUE)
      set(LIBINT_LIBRARIES ${LIBINT_LIBRARY})
      set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR})
      set(LIBINT_LIBRARY_DIRS ${LIBINT_LIBRARY_DIR})
      mark_as_advanced(LIBINT_INCLUDE_DIRS LIBINT_LIBRARY_DIRS LIBINT_LIBRARIES)
      message(STATUS "Found Libint libraries: ${LIBINT_LIBRARIES}")
    endif(LIBINT_LIBRARY)
  endif(LIBINT_INCLUDE_DIR)

message(STATUS "LIBINT INCLUDE DIRS FINAL = ${LIBINT_INCLUDE_DIRS}")

endif(LIBINT_INCLUDE_DIRS)
