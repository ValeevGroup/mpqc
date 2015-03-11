add_custom_target(examples)

set(sources
    sparse_integral_example.cpp
    sqrt_inverse_example.cpp
    hartree_fock.cpp
    lr_hartree_fock.cpp
)

# add_executable(sparse_integral_example EXCLUDE_FROM_ALL sparse_integral_example.cpp) 
# target_link_libraries(sparse_integral_example
#  ${TileClusterChem_LIBRARIES}
# )
# add_dependencies(examples sparse_integral_example)

add_executable(sqrt_inverse_example EXCLUDE_FROM_ALL sqrt_inverse_example.cpp) 
target_link_libraries(sqrt_inverse_example
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples sqrt_inverse_example)

add_executable(hartree_fock EXCLUDE_FROM_ALL hartree_fock.cpp) 
target_link_libraries(hartree_fock
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples hartree_fock)

add_executable(lr_hartree_fock EXCLUDE_FROM_ALL lr_hartree_fock.cpp) 
target_link_libraries(lr_hartree_fock
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples lr_hartree_fock)
