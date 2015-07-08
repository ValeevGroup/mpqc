add_custom_target(examples)
add_custom_target(sqrt_inverse)
add_custom_target(paper)
add_custom_target(bugs)

file(COPY ${CMAKE_SOURCE_DIR}/default_input.json DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/generate_input.py DESTINATION ${CMAKE_BINARY_DIR}/examples)

add_executable(hf_l_inv EXCLUDE_FROM_ALL hf_l_inv.cpp)
target_link_libraries(hf_l_inv
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(paper hf_l_inv)


add_executable(sqrt_inverse_example EXCLUDE_FROM_ALL sqrt_inverse_example.cpp) 
target_link_libraries(sqrt_inverse_example
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(sqrt_inverse sqrt_inverse_example)

add_executable(lr_hartree_fock EXCLUDE_FROM_ALL lr_hartree_fock.cpp) 
target_link_libraries(lr_hartree_fock
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples lr_hartree_fock)

add_executable(lr_hartree_fock_inv_oh EXCLUDE_FROM_ALL lr_hartree_fock_inv_oh.cpp) 
target_link_libraries(lr_hartree_fock_inv_oh
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples lr_hartree_fock_inv_oh)


add_executable(matrix_plot EXCLUDE_FROM_ALL matrix_plot.cpp) 
target_link_libraries(matrix_plot
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples matrix_plot)

add_executable(eigen_value_estimator_bug_min_example EXCLUDE_FROM_ALL 
    eigen_value_estimator_bug_min_example.cpp)
target_link_libraries(eigen_value_estimator_bug_min_example
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(bugs eigen_value_estimator_bug_min_example)
