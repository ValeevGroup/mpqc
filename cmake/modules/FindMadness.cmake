# - Try to find Madness lib
#
# Writen by Drew Lewis
# 

if (MADNESS_INCLUDE_DIRS)
    set(MADNESS_FOUND TRUE)

else()
    find_path(MADNESS_INCLUDE_DIR NAMES madness_config.h
         PATHS 
         ${MADNESS_DIR}/include 
         ${TILEDARRAY_DIR}/include
     )
       
    mark_as_advanced(MADNESS_INCLUDE_DIR)
    set(MADNESS_INCLUDE_DIRS ${MADNESS_INCLUDE_DIR} 
        CACHE PATH "The madness include path." )
        
    if(MADNESS_INCLUDE_DIRS)
        set(MADNESS_FOUND TRUE)
    endif()
  

endif(MADNESS_INCLUDE_DIRS)