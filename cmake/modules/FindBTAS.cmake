# Script to find the BTAS headers. 
set (BTAS_FOUND "NO")

if (BTAS_INSTALL_DIR)
    set (_BTAS_INSTALL_DIR ${BTAS_INSTALL_DIR})
endif(BTAS_INSTALL_DIR)

if (NOT _BTAS_INSTALL_DIR)
    message("ERROR: Unable able to find BTAS install directory. ${_BTAS_INSTALL_DIR}")
endif (NOT _BTAS_INSTALL_DIR)

set (BTAS_INC_SEARCH_DIR ${_BTAS_INSTALL_DIR})
find_path(BTAS_INCLUDE_DIR 
    btas/btas.h
    PATHS ${BTAS_INC_SEARCH_DIR} ENV CPATH
)
mark_as_advanced(BTAS_INCLUDE_DIR)

if (NOT BTAS_INCLUDE_DIR)
    message("ERROR: Unable to find BTAS include directory. Looked in ${BTAS_INC_SEARCH_DIR}.")
endif (NOT BTAS_INCLUDE_DIR)

if (BTAS_INCLUDE_DIR)
    set (BTAS_FOUND "YES")
    set (BTAS_INCLUDE_DIRS ${BTAS_INCLUDE_DIR} CACHE PATH "BTAS include directory")
    mark_as_advanced(BTAS_INCLUDE_DIRS)
endif(BTAS_INCLUDE_DIR)

if (BTAS_FOUND)
    message(STATUS "Found BTAS.")
endif(BTAS_FOUND)
