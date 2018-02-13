# -*- mode: cmake -*-

# Import Eigen provided by TiledArray

# make sure TiledAerray has been found
if (NOT TARGET TiledArray_Eigen)
  message (FATAL_ERROR "target TiledArray_Eigen not found, either find_package(TiledArray ...) has not been called yet or the discovered TiledArray is broken")
endif (NOT TARGET TiledArray_Eigen)

# import TiledArray_Eigen as MPQC_Eigen
add_library(MPQC_Eigen INTERFACE)
foreach(prop INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_LINK_LIBRARIES INTERFACE_POSITION_INDEPENDENT_CODE)
  get_property(EIGEN3_${prop} TARGET TiledArray_Eigen PROPERTY ${prop})
  set_property(TARGET MPQC_Eigen PROPERTY
          ${prop} ${EIGEN3_${prop}})
endforeach()
install(TARGETS MPQC_Eigen EXPORT mpqc COMPONENT mpqc)
