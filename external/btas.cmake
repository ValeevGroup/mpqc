# -*- mode: cmake -*-

# Import BTAS provided by TiledArray

# make sure TiledAerray has been found
if (NOT TARGET TiledArray_BTAS)
  message (FATAL_ERROR "target TiledArray_BTAS not found, either find_package(TiledArray ...) has not been called yet or the discovered TiledArray is broken")
endif (NOT TARGET TiledArray_BTAS)

# import TiledArray_BTAS as MPQC_BTAS
add_library(MPQC_BTAS INTERFACE)
foreach(prop INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_LINK_LIBRARIES INTERFACE_POSITION_INDEPENDENT_CODE)
  get_property(BTAS_${prop} TARGET TiledArray_BTAS PROPERTY ${prop})
  set_property(TARGET MPQC_BTAS PROPERTY
          ${prop} ${BTAS_${prop}})
endforeach()
install(TARGETS MPQC_BTAS EXPORT mpqc COMPONENT mpqc)
