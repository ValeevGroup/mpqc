add_custom_target(examples)
add_custom_target(fastK)

add_executable(cluster_example EXCLUDE_FROM_ALL cluster_example.cpp)
target_link_libraries(cluster_example molecule)
add_dependencies(examples cluster_example)

add_executable(atom_based_cluster_example EXCLUDE_FROM_ALL atom_based_cluster_example.cpp)
target_link_libraries(atom_based_cluster_example molecule)
add_dependencies(examples atom_based_cluster_example)

add_executable(molecule_example EXCLUDE_FROM_ALL molecule_example.cpp)
target_link_libraries(molecule_example molecule)
add_dependencies(examples molecule_example)

add_executable(clustering_example EXCLUDE_FROM_ALL clustering_example.cpp)
target_link_libraries(clustering_example molecule)
add_dependencies(examples clustering_example)

add_executable(basis_example EXCLUDE_FROM_ALL basis_example.cpp)
target_link_libraries(basis_example molecule;basis)
add_dependencies(examples basis_example)


add_executable(ccsd EXCLUDE_FROM_ALL ccsd.cpp)
target_link_libraries(ccsd ${TileClusterChem_LIBRARIES})

add_executable(mpqc EXCLUDE_FROM_ALL mpqc.cpp)
target_link_libraries(mpqc ${TileClusterChem_LIBRARIES})

add_executable(eomcc EXCLUDE_FROM_ALL eom_cc.cpp)
target_link_libraries(eomcc ${TileClusterChem_LIBRARIES})

file(COPY ${CMAKE_SOURCE_DIR}/examples/default_input.json DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/generate_input.py DESTINATION ${CMAKE_BINARY_DIR}/examples)

add_executable(just_world EXCLUDE_FROM_ALL empty_main.cpp) 
target_link_libraries(just_world
  ${TileClusterChem_LIBRARIES}
)

add_executable(new_ints EXCLUDE_FROM_ALL new_integrals.cpp) 
target_link_libraries(new_ints ${TileClusterChem_LIBRARIES})

add_executable(simple_direct_scf EXCLUDE_FROM_ALL simple_direct_scf.cpp) 
target_link_libraries(simple_direct_scf ${TileClusterChem_LIBRARIES})

add_executable(simple_stored_scf EXCLUDE_FROM_ALL simple_stored_scf.cpp) 
target_link_libraries(simple_stored_scf ${TileClusterChem_LIBRARIES})

add_executable(cadf_df_scf EXCLUDE_FROM_ALL cadf_df_scf.cpp) 
target_link_libraries(cadf_df_scf ${TileClusterChem_LIBRARIES})

add_executable(simple_stored_df_scf EXCLUDE_FROM_ALL simple_stored_df_scf.cpp) 
target_link_libraries(simple_stored_df_scf ${TileClusterChem_LIBRARIES})

add_executable(simple_direct_df_scf EXCLUDE_FROM_ALL simple_direct_df_scf.cpp) 
target_link_libraries(simple_direct_df_scf ${TileClusterChem_LIBRARIES})

add_executable(mp2_f12_expression EXCLUDE_FROM_ALL mp2_f12_expression.cpp)
target_link_libraries(mp2_f12_expression ${TileClusterChem_LIBRARIES})
