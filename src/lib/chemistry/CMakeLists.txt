set(subdirectories
  molecule
  qc
  solvent
)

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects ${${s}/OBJECTS})
endforeach()

list(APPEND objects $<TARGET_OBJECTS:molecule>)
list(APPEND objects $<TARGET_OBJECTS:solvent>)

add_library(chemistry ${objects})
set(chemistry/OBJECTS ${objects} PARENT_SCOPE)
install(TARGETS chemistry DESTINATION lib)

set(chemistry_deplibs math util)
if (HAVE_LIBINT2)
  set(chemistry_deplibs ${chemistry_deplibs} ${LIBINT2_LIBRARY})
endif()
if (HAVE_MADNESS)
  set(chemistry_deplibs ${chemistry_deplibs} ${MADNESS_LIBRARIES})
endif()
if (HAVE_PSI3)
  set(chemistry_deplibs ${chemistry_deplibs} ${PSI3_LIBRARY})
endif()
if (ARMCI_LIBRARIES)
  set(chemistry_deplibs ${chemistry_deplibs} ${ARMCI_LIBRARIES})
endif()
if (HAVE_BOOST)
  set(chemistry_deplibs ${chemistry_deplibs} ${Boost_LIBRARIES})
endif()
if (HAVE_HDF5)
  set(chemistry_deplibs ${chemistry_deplibs} ${HDF5_LIBRARIES})
endif()

target_link_libraries(chemistry ${chemistry_deplibs})
