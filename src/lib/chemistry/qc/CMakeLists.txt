set(subdirectories
  basis
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

if (MPQC_NEW_FEATURES_CI)
  list(APPEND subdirectories ci)
endif()

if (HAVE_LIBINT2)
  list(APPEND subdirectories libint2)
endif()

if (MPQC_R12)
  list(APPEND subdirectories mbptr12 ccr12)
endif()

if (HAVE_PSI3)
  list(APPEND subdirectories psi)
endif()

foreach(s ${subdirectories})
  add_subdirectory(${s})
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

set(qc/OBJECTS ${objects} PARENT_SCOPE)
