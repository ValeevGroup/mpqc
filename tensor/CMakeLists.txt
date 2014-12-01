set(SOURCES
    full_rank_tile.h
    tile_algebra.h
    tile_variant.h
    tile_ops.h
    tile_pimple.h
    low_rank_tile.h
    lr_tile.h
)

add_executable(flatten_test EXCLUDE_FROM_ALL flatten_test.cpp)
target_link_libraries(flatten_test ${MATH_LIBS})

add_subdirectory(examples)
add_subdirectory(eigen_tests)
add_subdirectory(profiling)
