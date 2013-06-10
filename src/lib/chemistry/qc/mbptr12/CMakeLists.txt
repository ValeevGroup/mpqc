# mbptr12 library

add_library(mbptr12 OBJECT
  FxFgen.cc
  Xgen.cc
  approxApp.cc
  approxB.cc
  approxC.cc
  approxC_GenRefansatz2.cc
  approxC_ansatz1.cc
  approxCp.cc
  bc_contrib.cc
  compute_VXB_a.cc
  compute_VXB_a_vbsneqobs.cc
  compute_VXB_abs.cc
  compute_VXB_c_ansatz1.cc
  compute_VX_GenRefansatz2.cc
  compute_amps.cc
  compute_density_diag.cc
  compute_energy_a.cc
  compute_tbint_tensor.cc
  coulomb.cc
  creator.cc
  dkh.cc
  ebc_contribs.cc
  exchange.cc
  fixedcoefficient.cc
  fock.cc
  gaussianfit.cc
  gbc_contribs.cc
  intermeds_g12.cc
  mbptr12.cc
  mp2_pair_energies.cc
  mp2r12_energy.cc
  mp2r12_energy_compute.cc
  mp2r12_energy_diag.cc
  mp2r12_energy_diag2.cc
  mp2r12_energy_util.cc
  multipole_ints.cc
  pt2r12.cc
  pt2r12_spinorbital.cc
  r12_amps.cc
  r12int_eval.cc
  r12technology.cc
  r12wfnworld.cc
  ref.cc
  ri_basis.cc
  singles_emp2.cc
  twobody_intermeds.cc
  twobodygrid.cc
  twoparticlecontraction.cc
)

# tests

#set_property(SOURCE mbptr12test.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(mbptr12test EXCLUDE_FROM_ALL mbptr12test.cc $<TARGET_OBJECTS:mbptr12>)
#add_test(mbptr12test mbptr12test)

