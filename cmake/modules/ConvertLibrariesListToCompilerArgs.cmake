#
# converts a list of libraries (second argument, don't forget to enclose the list in quotes)
# into a list of command-line parameters to the compiler/linker. The main purpose is to
# handle Apple's frameworks.
#

# transform library list _libs into compiler args _args
macro(convert_libs_to_compargs _args _libs )

 # for every member of _libs
 foreach (_lib ${_libs})
 
   if (TARGET "${_lib}") # ${_lib} = IMPORTED TARGET?
     
     get_property(_is_imported TARGET "${_lib}" PROPERTY IMPORTED)
     get_property(_lib_plus_deps_list TARGET "${_lib}" PROPERTY LOCATION)
     if (${_is_imported})
       # grab its dependents
       get_property(_lib_dep_libs TARGET "${_lib}" PROPERTY INTERFACE_LINK_LIBRARIES)
       list(APPEND _lib_plus_deps_list ${_lib_dep_libs})
       #message (STATUS "${_lib} is an IMPORTED TARGET located at ${_location}, _lib_plus_deps_list=${_lib_plus_deps_list}")
     endif(${_is_imported}) # skip non-IMPORTED targets
     # and process with convert_libs_to_compargs
     convert_libs_to_compargs(${_args} "${_lib_plus_deps_list}")
     
   else(TARGET "${_lib}") # ${_lib} = library?
   
     set(_libpath _libpath-NOTFOUND) # DON'T REMOVE THIS LINE
     find_library(_libpath ${_lib})
     if (_libpath)
       set(_library "${_libpath}")
     else(_libpath)
       set(_library "${_lib}")
     endif()
     #message(STATUS "_library = ${_library}")
     
     # Apple framework libs
     if (APPLE)
       get_filename_component(_ext ${_library} EXT)
       if ("${_ext}" STREQUAL ".framework")
         get_filename_component(_name ${_library} NAME_WE)
         get_filename_component(_path ${_library} PATH)
         #message(STATUS "${_library} = ${_name} ${_path} ${_ext}")
         set(_library "-F ${_path} -framework ${_name}")
       endif()
       #message(STATUS "${_lib} => ${_library}")
     endif(APPLE)
  
     set(${_args} "${${_args}} ${_library}")
     #message (STATUS "NOW ${_args}=${${_args}}")

   endif(TARGET "${_lib}")
  
 endforeach()
 #message (STATUS "convert_libs_to_compargs(_args=${_args},_libs=${_libs}) produced in ${_args}=${${_args}}")
 
endmacro()

