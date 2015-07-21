include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    cereal
    PREFIX ${CMAKE_BINARY_DIR}/include/cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (CEREAL_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(cereal source_dir)
set(CEREAL_INCLUDE_DIR ${source_dir}/include PARENT_SCOPE CACHE Path "Path to include folder for Cereal")
