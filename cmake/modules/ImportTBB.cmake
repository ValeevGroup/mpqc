# based on https://github.com/justusc/FindTBB
macro (import_tbb)
  if(TBB_FOUND)
  
    get_filename_component(TBB_LIBRARIES_DIR "${TBB_LIBRARY_DEBUG}" DIRECTORY)
    message(STATUS "Found Intel TBB: headers in ${TBB_INCLUDE_DIRS} (interface version ${TBB_INTERFACE_VERSION}), libs in ${TBB_LIBRARIES_DIR}")
    
    add_library(tbb SHARED IMPORTED)
    #message(STATUS "TBB_INCLUDE_DIRS=${TBB_INCLUDE_DIRS}")
    #message(STATUS "TBB_LIBRARIES=${TBB_LIBRARIES}")
    set_target_properties(tbb PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES  "${TBB_INCLUDE_DIRS}"
          IMPORTED_LOCATION              "${TBB_LIBRARIES}")
    if(TBB_LIBRARIES_RELEASE AND TBB_LIBRARIES_DEBUG)
      set_target_properties(tbb PROPERTIES
          INTERFACE_COMPILE_DEFINITIONS "$<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:TBB_USE_DEBUG=1>"
          IMPORTED_LOCATION_DEBUG          "${TBB_LIBRARY_DEBUG}"
          IMPORTED_LOCATION_RELWITHDEBINFO "${TBB_LIBRARY_DEBUG}"
          IMPORTED_LOCATION_RELEASE        "${TBB_LIBRARY_RELEASE}"
          IMPORTED_LOCATION_MINSIZEREL     "${TBB_LIBRARY_RELEASE}"
          )
    elseif(TBB_LIBRARIES_RELEASE)
      set_target_properties(tbb PROPERTIES IMPORTED_LOCATION "${TBB_LIBRARY_RELEASE}")
    else()
      set_target_properties(tbb PROPERTIES
          INTERFACE_COMPILE_DEFINITIONS "TBB_USE_DEBUG=1"
          IMPORTED_LOCATION              "${TBB_LIBRARY_DEBUG}"
          )
    endif()
  endif(TBB_FOUND)
endmacro()
