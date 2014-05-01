# - Try to find Elemental
#
# Writen by Drew Lewis
# 

if (ELEMENTAL_INCLUDE_DIRS)
    set(ELEMENTAL_FOUND TRUE)

else()
    find_path(ELEMENTAL_INCLUDE_DIR NAMES elemental.hpp
        PATHS ${ELEMENTAL_DIR}/include)
    
    mark_as_advanced(ELEMENTAL_INCLUDE_DIR)
    set(ELEMENTAL_INCLUDE_DIRS ${ELEMENTAL_INCLUDE_DIR}
        CACHE PATH "The Elemental include path.")
    
    if(ELEMENTAL_INCLUDE_DIRS)
        set(ELEMENTAL_FOUND TRUE)
    endif()
  
endif(ELEMENTAL_INCLUDE_DIRS)
