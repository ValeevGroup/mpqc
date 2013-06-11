set(subdirectories
  molecule
  qc
  solvent
)

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects ${${s}/OBJECTS})
endforeach()

list(APPEND objects $<TARGET_OBJECTS:molecule>)
list(APPEND objects $<TARGET_OBJECTS:solvent>)

add_library(chemistry STATIC ${objects})
set(chemistry/OBJECTS ${objects} PARENT_SCOPE)
install(TARGETS chemistry ARCHIVE DESTINATION lib)
