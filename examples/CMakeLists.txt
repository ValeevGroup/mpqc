add_custom_target(examples)
add_custom_target(read_density)

set(sources
    sparse_integral_example.cpp
    sqrt_inverse_example.cpp
    hartree_fock.cpp
    lr_hartree_fock.cpp
    four_center_exchange.cpp
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

add_executable(four_center_exchange EXCLUDE_FROM_ALL four_center_exchange.cpp) 
target_link_libraries(four_center_exchange
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(examples four_center_exchange)

# We need to save the old cmake cxx_flags and append boost include only for this 
# file.
set(CMAKE_CXX_FLAGS_ ${CMAKE_CXX_FLAGS}) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/${USER_BOOST_INCLUDE}")
add_executable(eri3_with_density_from_disk EXCLUDE_FROM_ALL eri3_with_density_from_disk.cpp) 
target_link_libraries(eri3_with_density_from_disk 
  ${TileClusterChem_LIBRARIES}
  ${USER_BOOST_LIB}
)

add_executable(eri3X EXCLUDE_FROM_ALL eri3X.cpp) 
target_link_libraries(eri3X
  ${TileClusterChem_LIBRARIES}
)
add_dependencies(read_density eri3_with_density_from_disk eri3X)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_}")
