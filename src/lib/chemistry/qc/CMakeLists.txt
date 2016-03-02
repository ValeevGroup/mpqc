include(AddLibrary)
SUBDIRLIST(Subdirs CMAKE_CURRENT_LIST_DIR)
FOREACH(subdir ${Subdirs})
    ADD_SUBDIRECTORY(${subdir})
ENDFOREACH()

set(subobjects
  basis
  dft
  etrain
  extern
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
  list(APPEND subobjects ci)
endif()

if (HAVE_LIBINT2)
  list(APPEND subobjects libint2)
endif()

if (MPQC_R12)
  list(APPEND subobjects mbptr12 ccr12)
endif()

if (HAVE_PSI3)
  list(APPEND subobjects psi)
endif()

add_mpqc_objects_target(qc ${subobjects})
