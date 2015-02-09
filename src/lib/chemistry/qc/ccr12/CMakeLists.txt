include(AddLibrary)
  
if (MPQC_R12)
  # ccr12 library
  # collect the .cc and .cpp files to be compiled.
  set(ccfiles
    ccr12.cc
    ccr12_info.cc
    ccr12_info_drivers.cc
    ccr12_info_drivers_r12.cc
    ccr12_info_fill_in.cc
    ccr12_info_fill_in_r12.cc
    ccr12_info_intermediates.cc
    ccr12_info_offset.cc
    ccr12_info_routines.cc
    ccr12_triples.cc
    ccsd.cc
    ccsd_2q_left.cc
    ccsd_2q_right.cc
    ccsd_2t_left.cc
    ccsd_2t_pr12_right.cc
    ccsd_2t_r12_left.cc
    ccsd_2t_right.cc
    ccsd_e.cc
    ccsd_pt.cc
    ccsd_pt_left.cc
    ccsd_pt_right.cc
    ccsd_r12.cc
    ccsd_r12_e.cc
    ccsd_r12_pt_right.cc
    ccsd_r12_t1.cc
    ccsd_r12_t2.cc
    ccsd_sub_bar_r12.cc
    ccsd_sub_r12_left.cc
    ccsd_sub_r12_right.cc
    ccsd_t1.cc
    ccsd_t2.cc
    ccsdpr12.cc
    ccsdpr12_c.cc
    ccsdpr12_t1.cc
    ccsdpr12_t2.cc
    ccsdt.cc
    ccsdt_t1.cc
    ccsdt_t2.cc
    ccsdt_t3.cc
    ccsdtq.cc
    ccsdtq_t2.cc
    ccsdtq_t3.cc
    ccsdtq_t4.cc
    lambda_ccsd_t1.cc
    lambda_ccsd_t2.cc
    lambda_ccsdpr12_t1.cc
    lambda_ccsdpr12_t2.cc
    mtensor.cc
    parenthesis2q.cc
    parenthesis2t.cc
    prediagon.cc
    tensor.cc
    tensorextrap.cc
    triples_denom_contraction.cc
    triples_doubles.cc
    triples_energy.cc
    triples_singles.cc
    triples_utils.cc
  )
  
  add_mpqc_object_library(ccr12 ${ccfiles})
  
  # tests
  
  if (MPQC_UNITTEST)
    # no tests yet
  endif(MPQC_UNITTEST)

endif (MPQC_R12)