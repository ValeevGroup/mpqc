set(subdirectories
  basis
  ccr12
  ci
#  cints
  dft
  etrain
  intv3
  lcao
  #libint2
  lmp2
  mbpt
  #mbptr12
  nbody
  oint3
  #psi
  scf
  wfn
  )

if (HAVE_LIBINT2)
  list(APPEND subdirectories libint2)
endif()

foreach(s ${subdirectories})
  add_subdirectory(${s})
#  message("${s} $<TARGET_OBJECTS:${s}>")
  list(APPEND objects $<TARGET_OBJECTS:${s}>)
endforeach()

set(qc/OBJECTS ${objects} PARENT_SCOPE)
