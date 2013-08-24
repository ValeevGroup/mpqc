include(AddLibrary)

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

if (MPQC_CI)
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

add_mpqc_objects_target(qc ${subdirectories})
