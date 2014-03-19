include(AddLibrary)

IF(HAVE_ELEMENTAL)
    add_mpqc_object_library(elemental
        init.cpp
        grid.cpp
    )
    
    IF(MPQC_UNITTEST)
        add_executable(elemtest element_test.cpp $<TARGET_OBJECTS:class>) 
        TARGET_LINK_LIBRARIES(elemtest util)
        TARGET_LINK_LIBRARIES(elemtest ${LAPACK_LIBRARIES})
        add_test(elemtest elemtest)
    ENDIF(MPQC_UNITTEST)
ENDIF(HAVE_ELEMENTAL)



