
add_subdirectory(lib)

add_dependencies(chemistry External)
add_dependencies(math External)
add_dependencies(util External)

add_subdirectory(bin)
