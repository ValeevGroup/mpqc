set(SOURCES
    full_rank_tile.h
    tile_algebra.h
    tile_variant.h
    tile_ops.h
    tile_pimple.h
    low_rank_tile.h
    lr_tile.h
)

add_subdirectory(examples)
add_subdirectory(eigen_tests)
add_subdirectory(profiling)
