//
// comp_g12nc.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdarg.h>
#include <cassert>

#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/libint2/macros.h>
#include <chemistry/qc/libint2/g12nc.h>
#include <chemistry/qc/libint2/tform.h>

using namespace std;
using namespace sc;

static inline void
swtch(GaussianBasisSet* &i,GaussianBasisSet* &j)
{
  GaussianBasisSet *tmp;
  tmp = i;
  i = j;
  j = tmp;
}

static inline void
pswtch(void**i,void**j)
{
  void*tmp;
  tmp = *i;
  *i = *j;
  *j = tmp;
}

static inline void
iswtch(int *i,int *j)
{
  int tmp;
  tmp = *i;
  *i = *j;
  *j = tmp;
}

static void
fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

void
G12NCLibint2::compute_quartet(int *psh1, int *psh2, int *psh3, int *psh4)
{
#ifdef EREP_TIMING
  char section[30];
#endif
  GaussianBasisSet *pbs1=bs1_.pointer();
  GaussianBasisSet *pbs2=bs2_.pointer();
  GaussianBasisSet *pbs3=bs3_.pointer();
  GaussianBasisSet *pbs4=bs4_.pointer();
  int int_expweight1; // For exponent weighted contractions.
  int int_expweight2;// For exponent weighted contractions.
  int int_expweight3;// For exponent weighted contractions.
  int int_expweight4;// For exponent weighted contractions.
  int size;
  int ii;
  int size1, size2, size3, size4;
  int tam1,tam2,tam3,tam4;
  int i,j,k,l;
  int pi, pj, pk, pl;
  int gci, gcj, gck, gcl;
  int sh1,sh2,sh3,sh4;
  int osh1,osh2,osh3,osh4;
  int am1,am2,am3,am4,am12,am34;
  int minam1,minam2,minam3,minam4;
  int redundant_index;
  int e12,e13e24,e34;
  int p12,p34,p13p24;
  int eAB;

#ifdef DMALLOC
  /*--- Test heap before ---*/
  int
  heapstate = dmalloc_verify(target_ints_buffer_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(cart_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(sphharm_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(perm_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(tformbuf_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
#endif

  osh1 = sh1 = *psh1;
  osh2 = sh2 = *psh2;
  osh3 = sh3 = *psh3;
  osh4 = sh4 = *psh4;

  /* Test the arguments to make sure that they are sensible. */
  if ( sh1 < 0 || sh1 >= bs1_->nbasis()
      || sh2 < 0 || sh2 >= bs2_->nbasis()
      || sh3 < 0 || sh3 >= bs3_->nbasis()
      || sh4 < 0 || sh4 >= bs4_->nbasis() ) {
    ExEnv::errn() << scprintf("compute_erep has been incorrectly used\n");
    ExEnv::errn() << scprintf("shells (bounds): %d (%d), %d (%d), %d (%d), %d (%d)\n",
        sh1,bs1_->nbasis()-1,
        sh2,bs2_->nbasis()-1,
        sh3,bs3_->nbasis()-1,
        sh4,bs4_->nbasis()-1);
    fail();
  }

  /* Set up pointers to the current shells. */
  int_shell1_ = &bs1_->shell(sh1);
  int_shell2_ = &bs2_->shell(sh2);
  int_shell3_ = &bs3_->shell(sh3);
  int_shell4_ = &bs4_->shell(sh4);

  /* Compute the maximum angular momentum on each centers to
   * determine the most efficient way to invoke the building and shifting
   * routines.  The minimum angular momentum will be computed at the
   * same time. */
  minam1 = int_shell1_->min_am();
  minam2 = int_shell2_->min_am();
  minam3 = int_shell3_->min_am();
  minam4 = int_shell4_->min_am();
  am1 = int_shell1_->max_am();
  am2 = int_shell2_->max_am();
  am3 = int_shell3_->max_am();
  am4 = int_shell4_->max_am();
  am12 = am1 + am2;
  am34 = am3 + am4;

  // This condition being true is guaranteed by the constructor of IntegralLibint2
  //if (minam1 != am1 ||
  //    minam2 != am2 ||
  //    minam3 != am3 ||
  //    minam4 != am4 ) {
  //  ExEnv::errn() << scprintf("Int2eLibint2::comp_g12() cannot yet handle fully general contractions") << endl;
  //  fail();
  //}
  
  /* See if need to transform to spherical harmonics */
  bool need_cart2sph_transform = false;
  if (int_shell1_->has_pure() ||
      int_shell2_->has_pure() ||
      int_shell3_->has_pure() ||
      int_shell4_->has_pure())
  need_cart2sph_transform = true;

  /* See if contraction quartets need to be resorted into a shell quartet */
  bool need_sort_to_shell_quartet = false;
  int num_gen_shells = 0;
  if (int_shell1_->ncontraction() > 1)
  num_gen_shells++;
  if (int_shell2_->ncontraction() > 1)
  num_gen_shells++;
  if (int_shell3_->ncontraction() > 1)
  num_gen_shells++;
  if (int_shell4_->ncontraction() > 1)
  num_gen_shells++;
  if (am12+am34 && num_gen_shells >= 1)
  need_sort_to_shell_quartet = true;

  /* Unique integrals are needed only ?*/
  bool need_unique_ints_only = false;
  if (!redundant_) {
    e12 = 0;
    if (int_shell1_ == int_shell2_ && int_shell1_->nfunction()>1)
    e12 = 1;
    e34 = 0;
    if (int_shell3_ == int_shell4_ && int_shell3_->nfunction()>1)
    e34 = 1;
    e13e24 = 0;
    if (int_shell1_ == int_shell3_ && int_shell2_ == int_shell4_ && int_shell1_->nfunction()*int_shell2_->nfunction()>1)
    e13e24 = 1;

    if ( e12 || e34 || e13e24 )
    need_unique_ints_only = true;
  }

#ifdef EREP_TIMING
  sprintf(section,"erep am=%02d",am12+am34);
  tim_enter(section);
  tim_enter("setup");
#endif

  /* Convert the integral to the most efficient form. */
  p12 = 0;
  p34 = 0;
  p13p24 = 0;

  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);
    iswtch(&minam1,&minam2);
    pswtch((void**)&int_shell1_,(void**)&int_shell2_);
    swtch(pbs1,pbs2);
  }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);
    iswtch(&minam3,&minam4);
    pswtch((void**)&int_shell3_,(void**)&int_shell4_);
    swtch(pbs3,pbs4);
  }
  if (am12 > am34) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);
    iswtch(&am12,&am34);
    iswtch(&minam1,&minam3);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell1_,(void**)&int_shell3_);
    swtch(pbs1,pbs3);
    pswtch((void**)&int_shell2_,(void**)&int_shell4_);
    swtch(pbs2,pbs4);
  }
  bool shells_were_permuted = (p12||p34||p13p24);

  /* If the centers were permuted, then the int_expweighted variable may
   * need to be changed. */
  if (p12) {
    iswtch(&int_expweight1,&int_expweight2);
  }
  if (p34) {
    iswtch(&int_expweight3,&int_expweight4);
  }
  if (p13p24) {
    iswtch(&int_expweight1,&int_expweight3);
    iswtch(&int_expweight2,&int_expweight4);
  }

  /* Compute the shell sizes. */
  size1 = int_shell1_->ncartesian();
  size2 = int_shell2_->ncartesian();
  size3 = int_shell3_->ncartesian();
  size4 = int_shell4_->ncartesian();
  size = size1*size2*size3*size4;

  /* Compute center data for Libint */
  int ctr1 = pbs1->shell_to_center(sh1);
  int ctr2 = pbs2->shell_to_center(sh2);
  int ctr3 = pbs3->shell_to_center(sh3);
  int ctr4 = pbs4->shell_to_center(sh4);
  Libint_t& libint0 = Libint_[0];
  libint0.AB_x[0] = pbs1->r(ctr1,0) - pbs2->r(ctr2,0);
  libint0.AB_y[0] = pbs1->r(ctr1,1) - pbs2->r(ctr2,1);
  libint0.AB_z[0] = pbs1->r(ctr1,2) - pbs2->r(ctr2,2);
  libint0.CD_x[0] = pbs3->r(ctr3,0) - pbs4->r(ctr4,0);
  libint0.CD_y[0] = pbs3->r(ctr3,1) - pbs4->r(ctr4,1);
  libint0.CD_z[0] = pbs3->r(ctr3,2) - pbs4->r(ctr4,2);
  for(i=0;i<3;i++) {
    quartet_info_.A[i] = pbs1->r(ctr1,i);
    quartet_info_.B[i] = pbs2->r(ctr2,i);
    quartet_info_.C[i] = pbs3->r(ctr3,i);
    quartet_info_.D[i] = pbs4->r(ctr4,i);
  }
  quartet_info_.AB2 = libint0.AB_x[0]*libint0.AB_x[0];
  quartet_info_.AB2 += libint0.AB_y[0]*libint0.AB_y[0];
  quartet_info_.AB2 += libint0.AB_z[0]*libint0.AB_z[0];
  quartet_info_.CD2 = libint0.CD_x[0]*libint0.CD_x[0];
  quartet_info_.CD2 += libint0.CD_y[0]*libint0.CD_y[0];
  quartet_info_.CD2 += libint0.CD_z[0]*libint0.CD_z[0];

  /* Set up pointers to the current shell pairs. */
  quartet_info_.shell_pair12 = shell_pairs12_->shell_pair(osh1,osh2);
  quartet_info_.shell_pair34 = shell_pairs34_->shell_pair(osh3,osh4);

  /* Remember how permuted - will need to access shell pairs in grt_quartet_data_() using the original
   primitive indices */
  quartet_info_.p12 = p12;
  quartet_info_.p34 = p34;
  quartet_info_.p13p24 = p13p24;

  /* Remember the original primitive indices to access shell pair data
   Note the reverse order of switching, p13p24 first,
   then p12 and p34 - because we need the inverse mapping! */
  quartet_info_.op1 = &quartet_info_.p1;
  quartet_info_.op2 = &quartet_info_.p2;
  quartet_info_.op3 = &quartet_info_.p3;
  quartet_info_.op4 = &quartet_info_.p4;
  if (p13p24) {
    pswtch((void **)&quartet_info_.op1,(void **)&quartet_info_.op3);
    pswtch((void **)&quartet_info_.op2,(void **)&quartet_info_.op4);
  }
  if (p12)
  pswtch((void **)&quartet_info_.op1,(void **)&quartet_info_.op2);
  if (p34)
  pswtch((void **)&quartet_info_.op3,(void **)&quartet_info_.op4);

  /* Determine where integrals need to go at each stage */
  if (shells_were_permuted)
  if (need_sort_to_shell_quartet) {
    for(int te_type=0; te_type<num_te_types_; te_type++)
    prim_ints_[te_type] = cart_ints_[te_type];
    if (need_cart2sph_transform)
    for(int te_type=0; te_type<num_te_types_; te_type++)
    contr_quartets_[te_type] = sphharm_ints_;
    else
    for(int te_type=0; te_type<num_te_types_; te_type++)
    contr_quartets_[te_type] = cart_ints_[te_type];
    for(int te_type=0; te_type<num_te_types_; te_type++)
    shell_quartet_[te_type] = perm_ints_;
  }
  else {
    for(int te_type=0; te_type<num_te_types_; te_type++)
    prim_ints_[te_type] = cart_ints_[te_type];
    if (need_cart2sph_transform) {
      for(int te_type=0; te_type<num_te_types_; te_type++) {
        contr_quartets_[te_type] = sphharm_ints_;
        shell_quartet_[te_type] = contr_quartets_[te_type];
      }
    }
    else
    for(int te_type=0; te_type<num_te_types_; te_type++)
    shell_quartet_[te_type] = cart_ints_[te_type];
  }
  else
  if (need_sort_to_shell_quartet) {
    for(int te_type=0; te_type<num_te_types_; te_type++)
    prim_ints_[te_type] = cart_ints_[te_type];
    if (need_cart2sph_transform)
    for(int te_type=0; te_type<num_te_types_; te_type++)
    contr_quartets_[te_type] = sphharm_ints_;
    else
    for(int te_type=0; te_type<num_te_types_; te_type++)
    contr_quartets_[te_type] = cart_ints_[te_type];
    for(int te_type=0; te_type<num_te_types_; te_type++)
    shell_quartet_[te_type] = target_ints_buffer_[te_type];
  }
  else {
    if (need_cart2sph_transform) {
      for(int te_type=0; te_type<num_te_types_; te_type++) {
        prim_ints_[te_type] = cart_ints_[te_type];
        contr_quartets_[te_type] = target_ints_buffer_[te_type];
        shell_quartet_[te_type] = target_ints_buffer_[te_type];
      }
    }
    else {
      for(int te_type=0; te_type<num_te_types_; te_type++) {
        prim_ints_[te_type] = target_ints_buffer_[te_type];
        shell_quartet_[te_type] = target_ints_buffer_[te_type];
      }
    }
  }

  //
  // Need to distinguish 1-geminal and 2-geminal cases, hence don't use geminal_xxx_ directly, via references only
  //
  const bool braonly = (geminal_ket_ == IntParamsG12::null_geminal);
  const ContractedGeminal& gbra = geminal_bra_;
  const ContractedGeminal& gket = braonly ? IntParamsG12::zero_exponent_geminal : geminal_ket_;
  const int nbra_geminal_prims = gbra.size();
  const int nket_geminal_prims = gket.size();

  // Compute integrals one type at a time
  for(int te_type = 0; te_type < 4; te_type++) {

    OperType otype;
    switch (te_type) {
      case 0: otype = coulomb; break;
      case 1: otype = braonly ? f12 : f12_2; break;
      case 2: otype = f12_coulomb; break;
      case 3: otype = f12_T1_f12; break;
      default: MPQC_ASSERT(false); // unreachable
    }
    // no need f12/r12 when have 2 geminals
    if (otype == f12_coulomb && !braonly)
      continue;
    // need f12,T1,f12 only when have 2 geminals
    if (otype == f12_T1_f12 && braonly)
      continue;

    /* Begin loops over generalized contractions. */
    int buffer_offset = 0;
    for (gci=0; gci<int_shell1_->ncontraction(); gci++) {
      tam1 = int_shell1_->am(gci);
      int tsize1 = INT_NCART_NN(tam1);
      quartet_info_.gc1 = gci;
      for (gcj=0; gcj<int_shell2_->ncontraction(); gcj++) {
        tam2 = int_shell2_->am(gcj);
        int tsize2 = INT_NCART_NN(tam2);
        quartet_info_.gc2 = gcj;
        for (gck=0; gck<int_shell3_->ncontraction(); gck++) {
          tam3 = int_shell3_->am(gck);
          int tsize3 = INT_NCART_NN(tam3);
          quartet_info_.gc3 = gck;
          for (gcl=0; gcl<int_shell4_->ncontraction(); gcl++) {
            tam4 = int_shell4_->am(gcl);
            int tsize4 = INT_NCART_NN(tam4);
            quartet_info_.gc4 = gcl;
            quartet_info_.am = tam1 + tam2 + tam3 + tam4;
            const int size = tsize1*tsize2*tsize3*tsize4;

            // zero out the contracted integrals buffer
            memset(&(prim_ints_[te_type][buffer_offset]),0,size*sizeof(double));
#  if LIBINT2_ACCUM_INTS
            // zero out targets in Libint_
            Libint_[0].zero_out_targets = 1;
#  endif

            // Begin loop over basis function primitives
            int num_prim_combinations = 0;
            for (pi=0; pi<int_shell1_->nprimitive(); pi++) {
              quartet_info_.p1 = pi;
              for (pj=0; pj<int_shell2_->nprimitive(); pj++) {
                quartet_info_.p2 = pj;
                for (pk=0; pk<int_shell3_->nprimitive(); pk++) {
                  quartet_info_.p3 = pk;
                  for (pl=0; pl<int_shell4_->nprimitive(); pl++) {
                    quartet_info_.p4 = pl;

                    g12nc_quartet_data_(&(Libint_[num_prim_combinations]), 1.0, otype,
                                        &gbra, &gket);

                    if (quartet_info_.am == 0) {
                      prim_ints_[te_type][buffer_offset] += Libint_[num_prim_combinations].LIBINT_T_SS_EREP_SS(0)[0];
                    }

                    ++num_prim_combinations;
                  }
                }
              }
            }

            // compute geminal integrals
            if (quartet_info_.am) {
              // Compute the integrals
              Libint_[0].contrdepth = num_prim_combinations;
              LIBINT2_PREFIXED_NAME(libint2_build_eri)[tam1][tam2][tam3][tam4](&Libint_[0]);
              const LIBINT2_REALTYPE* prim_ints = Libint_[0].targets[0];
              std::copy(prim_ints, prim_ints+size, prim_ints_[te_type]+buffer_offset);
            }
            else {
              // ssss integrals were computed in the primitive loop
            }

            buffer_offset += size;
          }}}} // loops over contractions
    } // loop over operator types

  for(int te_type=0; te_type < num_te_types_; te_type++) {
    /*-------------------------------------------
     Transform to spherical harmonics if needed
     -------------------------------------------*/
    if (need_cart2sph_transform)
    transform_contrquartets_(prim_ints_[te_type],contr_quartets_[te_type]);

    //
    // If not CCA-compliant normalization -- re-normalize all integrals to 1
    //
#if INTEGRALLIBINT2_NORMCONV != INTEGRALLIBINT2_NORMCONV_CCA
    norm_contrcart1_(need_cart2sph_transform ? contr_quartets_[te_type] : prim_ints_[te_type]);
#endif

    /*----------------------------------------------
     Resort integrals from by-contraction-quartets
     into shell-quartet order if needed
     ----------------------------------------------*/
    if (need_sort_to_shell_quartet)
    sort_contrquartets_to_shellquartet_(contr_quartets_[te_type],shell_quartet_[te_type]);

    /*---------------------------------
     Permute integrals back if needed
     ---------------------------------*/
    if ((!permute_)&&shells_were_permuted) {
      // handle integrals first
      permute_target_(shell_quartet_[te_type],target_ints_buffer_[te_type],p13p24,p12,p34);
    }
  }
  if ((!permute_)&&shells_were_permuted) {
    // then indices
    if (p13p24) {
      iswtch(&sh1,&sh3);iswtch(psh1,psh3);
      iswtch(&sh2,&sh4);iswtch(psh2,psh4);
      iswtch(&am1,&am3);
      iswtch(&am2,&am4);
      iswtch(&am12,&am34);
      pswtch((void**)&int_shell1_,(void**)&int_shell3_);
      swtch(pbs1,pbs3);
      pswtch((void**)&int_shell2_,(void**)&int_shell4_);
      swtch(pbs2,pbs4);
      iswtch(&int_expweight1,&int_expweight3);
      iswtch(&int_expweight2,&int_expweight4);
    }
    if (p34) {
      iswtch(&sh3,&sh4);iswtch(psh3,psh4);
      iswtch(&am3,&am4);
      pswtch((void**)&int_shell3_,(void**)&int_shell4_);
      swtch(pbs3,pbs4);
      iswtch(&int_expweight3,&int_expweight4);
    }
    if (p12) {
      iswtch(&sh1,&sh2);iswtch(psh1,psh2);
      iswtch(&am1,&am2);
      pswtch((void**)&int_shell1_,(void**)&int_shell2_);
      swtch(pbs1,pbs2);
      iswtch(&int_expweight1,&int_expweight2);
    }
  }

  // Don't know how to extract nonredundant integrals because of non-Hermitian operators
  if (need_unique_ints_only) {
    throw FeatureNotImplemented("G12Libint2::compute_quartet() -- need_unique_ints_only is not implemented for this evaluator",__FILE__,__LINE__);
  }

#ifdef DMALLOC
  /*--- Test heap after ---*/
  heapstate = dmalloc_verify(target_ints_buffer_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(cart_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(sphharm_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(perm_ints_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
  heapstate = dmalloc_verify(tformbuf_);
  if (heapstate == DMALLOC_VERIFY_ERROR)
  fail();
#endif

  return;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
