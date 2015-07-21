include(ExternalProject)
find_package(Git REQUIRED)
add_custom_target(unit_tests)

ExternalProject_Add(
    catch
    PREFIX ${CMAKE_BINARY_DIR}/tests/catch
    GIT_REPOSITORY https://github.com/philsquared/Catch.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

set_target_properties(catch PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch source_dir)
set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE Path "Path to include folder for Catch")
include_directories(${CATCH_INCLUDE_DIR})

add_dependencies(unit_tests catch)

set(utests
    atom_test.cpp
    cluster_test.cpp
    )

add_executable(tests EXCLUDE_FROM_ALL test_main.cpp ${utests})
target_link_libraries(tests ${TileClusterChem_LIBRARIES})
set_target_properties(unit_tests PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_dependencies(tests catch)
add_dependencies(unit_tests tests)
add_test(NAME TestAll COMMAND ${CMAKE_BINARY_DIR}/tests/tests)
