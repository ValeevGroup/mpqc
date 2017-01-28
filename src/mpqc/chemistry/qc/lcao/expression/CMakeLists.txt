
set(sources
        operator.h
        operator.cpp
        orbital_index.h
        orbital_index.cpp
        formula.h
        formula.cpp
        formula_registry.h
        greek_to_english_name.h
        orbital_space.h
        orbital_registry.h
        permutation.h
        permutation.cpp
        )

add_mpqc_library(lcao_expression sources sources "MPQCmath_ta" "mpqc/chemistry/qc/lcao/expression")
