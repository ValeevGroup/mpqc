set(subdirectories
  basis
#  ci
  dft
  etrain
  intv3
  lcao
  lmp2
  mbpt
  nbody
  oint3
  scf
  wfn
  )

if (HAVE_LIBINT2)
  list(APPEND subdirectories libint2 mbptr12 ccr12)
endif()

if (HAVE_PSI3)
  list(APPEND subdirectories psi)
endif()

foreach(s ${subdirectories})
  add_subdirectory(${s})
#  message("${s} $<TARGET_OBJECTS:${s}>")
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

set(qc/OBJECTS ${objects} PARENT_SCOPE)
