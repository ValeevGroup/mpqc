# - Try to find TiledArray
#
# Writen by Drew Lewis
# 

if (TILEDARRAY_INCLUDE_DIRS)
    set(TILEDARRAY_FOUND TRUE)

else()
    find_path(TILEDARRAY_INCLUDE_DIR NAMES tiledarray.h 
        PATHS ${TILEDARRAY_DIR}/include)
    
    mark_as_advanced(TILEDARRAY_INCLUDE_DIR)
    
    set(TILEDARRAY_INCLUDE_DIRS ${TILEDARRAY_INCLUDE_DIR}
        CACHE PATH "The TiledArray include path.")

    if(TILEDARRAY_INCLUDE_DIRS)
        SET(TILEDARRAY_FOUND TRUE)
    endif()
  
endif(TILEDARRAY_INCLUDE_DIRS)