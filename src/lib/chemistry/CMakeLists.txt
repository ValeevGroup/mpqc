include(AddLibrary)

SUBDIRLIST(Subdirs CMAKE_CURRENT_SOURCE_DIR)

FOREACH(subdir ${Subdirs})
    ADD_SUBDIRECTORY(${subdir})
ENDFOREACH()

set(chemobjects
  molecule
  qc
  solvent
)

add_mpqc_library(chemistry ${chemobjects})
add_dependencies(math util)

target_link_libraries(chemistry math util)

if (HAVE_OPENBABEL2)
  target_link_libraries(chemistry ${OPENBABEL2_LIBRARIES})
endif()

if (HAVE_LIBINT2)
  target_link_libraries(chemistry ${LIBINT2_LIBRARY})
endif()

if (HAVE_PSI3)
  target_link_libraries(chemistry ${PSI3_LIBRARIES})
endif()

if (HAVE_BOOST)
  target_link_libraries(chemistry ${Boost_LIBRARIES})
  if (PYTHONLIBS_FOUND)
    target_link_libraries(chemistry ${PYTHON_LIBRARIES})
  endif()
endif()

if (HAVE_HDF5)
  target_link_libraries(chemistry ${HDF5_LIBRARIES})
endif()

if (HAVE_ARMCI)
  target_link_libraries(chemistry ${ARMCI_LIBRARIES})
endif()

if (HAVE_MPI)
  target_link_libraries(chemistry ${MPI_LIBRARIES})
endif()

install(TARGETS chemistry DESTINATION lib)
