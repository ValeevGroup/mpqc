
set(sources
        operation.h
        operation.cpp
        orbital_index.h
        orbital_index.cpp
        formula.h
        formula.cpp
        formula_registry.h
        greek_to_english_name.h
        orbital_space.h
        orbital_registry.h
        )

add_library(expression ${sources})
target_link_libraries(expression ${TA_LIBRARIES})
