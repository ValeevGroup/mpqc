add_custom_target(examples)
add_custom_target(sqrt_inverse)
add_custom_target(paper)
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

add_executable(mp2 EXCLUDE_FROM_ALL mp2.cpp)
target_link_libraries(mp2 ${TileClusterChem_LIBRARIES})

add_executable(ccsd EXCLUDE_FROM_ALL ccsd.cpp)
target_link_libraries(ccsd ${TileClusterChem_LIBRARIES})

file(COPY ${CMAKE_SOURCE_DIR}/examples/default_input.json DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/generate_input.py DESTINATION ${CMAKE_BINARY_DIR}/examples)

add_executable(hf_l_inv EXCLUDE_FROM_ALL hf_l_inv.cpp)
target_link_libraries(hf_l_inv
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(paper hf_l_inv)

add_executable(hf_l_inv_trick EXCLUDE_FROM_ALL hf_l_inv_head_gordon_trick.cpp)
target_link_libraries(hf_l_inv_trick
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(fastK hf_l_inv_trick)

add_executable(sqrt_inverse_example EXCLUDE_FROM_ALL sqrt_inverse_example.cpp) 
target_link_libraries(sqrt_inverse_example
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(sqrt_inverse sqrt_inverse_example)

add_executable(matrix_plot EXCLUDE_FROM_ALL matrix_plot.cpp) 
target_link_libraries(matrix_plot
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples matrix_plot)

add_executable(just_world EXCLUDE_FROM_ALL empty_main.cpp) 
target_link_libraries(just_world
  ${TileClusterChem_LIBRARIES}
)

add_executable(old_ints EXCLUDE_FROM_ALL old_integrals.cpp) 
target_link_libraries(old_ints
  ${TileClusterChem_LIBRARIES}
)

add_executable(new_ints EXCLUDE_FROM_ALL new_integrals.cpp) 
target_link_libraries(new_ints ${TileClusterChem_LIBRARIES})

