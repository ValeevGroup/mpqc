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
target_link_libraries(chemistry ${OPENBABEL2_LIBRARIES})
target_link_libraries(chemistry ${LIBINT2_LIBRARY})
target_link_libraries(chemistry ${PSI3_LIBRARIES})
target_link_libraries(chemistry ${Boost_LIBRARIES})
target_link_libraries(chemistry ${PYTHON_LIBRARIES})
target_link_libraries(chemistry ${HDF5_LIBRARIES})
target_link_libraries(chemistry ${ARMCI_LIBRARIES})
target_link_libraries(chemistry ${MPI_LIBRARIES})

install(TARGETS chemistry DESTINATION lib)
