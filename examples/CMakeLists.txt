add_custom_target(examples)
add_custom_target(paper)
add_custom_target(read_density)

add_executable(hf_for_paper EXCLUDE_FROM_ALL hf_for_paper.cpp)
target_link_libraries(hf_for_paper
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(paper hf_for_paper)

add_executable(hf_l_inv EXCLUDE_FROM_ALL hf_l_inv.cpp)
target_link_libraries(hf_l_inv
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(paper hf_l_inv)


add_executable(sqrt_inverse_example EXCLUDE_FROM_ALL sqrt_inverse_example.cpp) 
target_link_libraries(sqrt_inverse_example
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples sqrt_inverse_example)

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

add_executable(new_tile_test EXCLUDE_FROM_ALL new_tile_test.cpp)
target_link_libraries(new_tile_test
    ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples new_tile_test)

add_executable(matrix_plot EXCLUDE_FROM_ALL matrix_plot.cpp) 
target_link_libraries(matrix_plot
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples matrix_plot)

add_executable(eri3X EXCLUDE_FROM_ALL eri3X.cpp) 
target_link_libraries(eri3X
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples eri3X)

# We need to save the old cmake cxx_flags and append boost include only for this 
# file.
set(CMAKE_CXX_FLAGS_ ${CMAKE_CXX_FLAGS}) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/${USER_BOOST_INCLUDE}")
add_executable(eri3_with_density_from_disk EXCLUDE_FROM_ALL eri3_with_density_from_disk.cpp) 
target_link_libraries(eri3_with_density_from_disk 
  ${TileClusterChem_LIBRARIES}
  ${USER_BOOST_LIB}
)


add_executable(eri3X12 EXCLUDE_FROM_ALL eri3X12.cpp) 
target_link_libraries(eri3X12
  ${TileClusterChem_LIBRARIES}
  ${USER_BOOST_LIB}
)

add_executable(eri3X_with_density EXCLUDE_FROM_ALL eri3X_with_density.cpp) 
target_link_libraries(eri3X_with_density
  ${TileClusterChem_LIBRARIES}
  ${USER_BOOST_LIB}
)

add_dependencies(read_density eri3_with_density_from_disk eri3X12 eri3X_with_density)

# Fix cmake flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_}")
