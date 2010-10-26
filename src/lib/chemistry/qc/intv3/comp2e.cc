//
// comp2e.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#include <util/misc/formio.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/tformv3.h>

using namespace std;
using namespace sc;

#undef DER_TIMING
#undef EREP_TIMING

#if defined(DER_TIMING)||defined(EREP_TIMING)
#  include <util/misc/regtime.h>
#endif

static inline void
swtch(GaussianBasisSet* &i,GaussianBasisSet* &j)
{
  GaussianBasisSet *tmp;
  tmp = i;
  i = j;
  j = tmp;
}

static inline void
sswtch(GaussianShell**i,GaussianShell**j)
{
  GaussianShell*tmp;
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

/* This computes the 2erep integrals for a shell quartet
 * specified by psh1, psh2, psh3, psh4.
 * The routine int_initialize_2erep must be called before
 * any integrals can be computed.
 * This routine may decide to change the shell ordering.
 * The new ordering is placed in *psh1,4 on exit.
 * for the derivatives.
 */
void
Int2eV3::erep(int &psh1, int &psh2, int &psh3, int &psh4)
{
  compute_erep(0,&psh1,&psh2,&psh3,&psh4,0,0,0,0);
  }

/* This is an alternate interface to int_erep.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
void
Int2eV3::erep(int *shells, int  *sizes)
{
  erep(shells[0],shells[1],shells[2],shells[3]);
  if (sizes) {
    sizes[0] = bs1_->shell(shells[0]).nfunction();
    sizes[1] = bs2_->shell(shells[1]).nfunction();
    sizes[2] = bs3_->shell(shells[2]).nfunction();
    sizes[3] = bs4_->shell(shells[3]).nfunction();
    }
  }

/* If we need a computation with adjusted angular momentum, then
 * this lower level routine can be called instead of int_erep.
 * The dam{1,2,3,4} arguments given the amount by which the
 * angular momentum is to adjusted.  This differs from libint version
 * 1 in that it used total angular momentum here.
 */
void
Int2eV3::compute_erep(int flags, int *psh1, int *psh2, int *psh3, int *psh4,
                      int dam1, int dam2, int dam3, int dam4)
{
#ifdef EREP_TIMING
  char section[30];
#endif
  GaussianBasisSet *pbs1=bs1_.pointer();
  GaussianBasisSet *pbs2=bs2_.pointer();
  GaussianBasisSet *pbs3=bs3_.pointer();
  GaussianBasisSet *pbs4=bs4_.pointer();
  int size;
  int ii;
  int size1, size2, size3, size4;
  int tam1,tam2,tam3,tam4;
  int i,j,k,l;
  int ogc1,ogc2,ogc3,ogc4;
  int sh1,sh2,sh3,sh4;
  int am1,am2,am3,am4,am12, am34;
  int minam1,minam2,minam3,minam4;
  int redundant_index;
  int e12,e13e24,e34;
  int p12,p34,p13p24;
  int eAB;

  /* Compute the offset shell numbers. */
  osh1 = *psh1 + bs1_shell_offset_;
  osh2 = *psh2 + bs2_shell_offset_;
  osh3 = *psh3 + bs3_shell_offset_;
  osh4 = *psh4 + bs4_shell_offset_;

  sh1 = *psh1;
  sh2 = *psh2;
  sh3 = *psh3;
  sh4 = *psh4;

  /* Test the arguments to make sure that they are sensible. */
  if (   sh1 < 0 || sh1 >= bs1_->nbasis()
      ||( !int_unit2 && (sh2 < 0 || sh2 >= bs2_->nbasis()))
      || sh3 < 0 || sh3 >= bs3_->nbasis()
      ||( !int_unit4 && (sh4 < 0 || sh4 >= bs4_->nbasis()))) {
    ExEnv::errn() << scprintf("compute_erep has been incorrectly used\n");
    ExEnv::errn() << scprintf("shells (bounds): %d (%d), %d (%d), %d (%d), %d (%d)\n",
            sh1,bs1_->nbasis()-1,
            sh2,(bs2_.null()?0:bs2_->nbasis())-1,
            sh3,bs3_->nbasis()-1,
            sh4,(bs4_.null()?0:bs4_->nbasis())-1);
    fail();
    }

  /* Set up pointers to the current shells. */
  int_shell1 = &bs1_->shell(sh1);
  if (!int_unit2) int_shell2 = &bs2_->shell(sh2);
  else int_shell2 = int_unit_shell;
  int_shell3 = &bs3_->shell(sh3);
  if (!int_unit4) int_shell4 = &bs4_->shell(sh4);
  else int_shell4 = int_unit_shell;


  /* Compute the maximum angular momentum on each centers to
   * determine the most efficient way to invoke the building and shifting
   * routines.  The minimum angular momentum will be computed at the
   * same time. */
  minam1 = int_shell1->min_am();
  minam2 = int_shell2->min_am();
  minam3 = int_shell3->min_am();
  minam4 = int_shell4->min_am();
  am1 = int_shell1->max_am();
  am2 = int_shell2->max_am();
  am3 = int_shell3->max_am();
  am4 = int_shell4->max_am();

  am1 += dam1; minam1 += dam1;
  am2 += dam2; minam2 += dam2;
  am3 += dam3; minam3 += dam3;
  am4 += dam4; minam4 += dam4;
  am12 = am1 + am2;
  am34 = am3 + am4;

  /* if no angular momentum remains on one of the centers return */
  if (am1 < 0 || am2 < 0 || am3 < 0 || am4 < 0) {
    return;
    }

#ifdef EREP_TIMING
  sprintf(section,"erep am=%02d",am12+am34);
  Timer tim(section);
  tim.enter("setup");
#endif

  /* Convert the integral to the most efficient form. */
  p12 = 0;
  p34 = 0;
  p13p24 = 0;

  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
    iswtch(&dam1,&dam2);
    iswtch(&minam1,&minam2);
    sswtch(&int_shell1,&int_shell2);
    swtch(pbs1,pbs2);
    }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
    iswtch(&dam3,&dam4);
    iswtch(&minam3,&minam4);
    sswtch(&int_shell3,&int_shell4);
    swtch(pbs3,pbs4);
    }
  if ((osh1 == osh4) && (osh2 == osh3) && (osh1 != osh2)) {
    /* Don't make the permutation unless we won't override what was
     * decided above about p34. */
    if (am4 == am3) {
      p34 = 1;
      iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&dam3,&dam4);
      iswtch(&minam3,&minam4);
      sswtch(&int_shell3,&int_shell4);
      swtch(pbs3,pbs4);
      }
    }
  if ((am34 > am12)||((am34 == am12)&&(minam1 > minam3))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&int_unit2,&int_unit4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    sswtch(&int_shell1,&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    sswtch(&int_shell2,&int_shell4);
    swtch(pbs2,pbs4);
    }
  /* This tries to make centers A and B equivalent, if possible. */
  else if (  (am3 == am1)
           &&(am4 == am2)
           && !int_unit2
           && !int_unit4
           &&(minam1 == minam3)
           &&(!(  (bs1_ == bs2_)
                &&(bs1_->shell_to_center(sh1)==bs2_->shell_to_center(sh2))))
           &&(   (bs3_ == bs4_)
               &&(bs3_->shell_to_center(sh3)==bs4_->shell_to_center(sh4)))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    sswtch(&int_shell1,&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    sswtch(&int_shell2,&int_shell4);
    swtch(pbs2,pbs4);
    }

  if ((pbs1 == pbs2)
      &&(pbs1->shell_to_center(sh1)==pbs2->shell_to_center(sh2))) {
    eAB = 1;
    }
  else {
    eAB = 0;
    }

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

  pbs1_ = pbs1;
  pbs2_ = pbs2;
  pbs3_ = pbs3;
  pbs4_ = pbs4;

  int nc1 = int_shell1->ncontraction();
  int nc2 = int_shell2->ncontraction();
  int nc3 = int_shell3->ncontraction();
  int nc4 = int_shell4->ncontraction();

  /* Compute the shell sizes. */
  for (ii=size1=0; ii<nc1; ii++)
    size1 += INT_NCART(int_shell1->am(ii)+dam1);
  for (ii=size2=0; ii<nc2; ii++)
    size2 += INT_NCART(int_shell2->am(ii)+dam2);
  for (ii=size3=0; ii<nc3; ii++)
    size3 += INT_NCART(int_shell3->am(ii)+dam3);
  for (ii=size4=0; ii<nc4; ii++)
    size4 += INT_NCART(int_shell4->am(ii)+dam4);
  size = size1*size2*size3*size4;

  if (int_integral_storage) {
#ifdef EREP_TIMING
      tim.change("check storage");
#endif
    if (dam1 || dam2 || dam3 || dam4) {
      ExEnv::errn() << scprintf("cannot use integral storage and dam\n");
      fail();
      }
    if (    !int_unit2
         && !int_unit4
         && int_have_stored_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24)) {
      goto post_computation;
      }
    }

  /* Buildam up on center 1 and 3. */
#ifdef EREP_TIMING
  tim.change("build");
#endif
  int_buildgcam(minam1,minam2,minam3,minam4,
                am1,am2,am3,am4,
                dam1,dam2,dam3,dam4,
                sh1,sh2,sh3,sh4,
                eAB);
#ifdef EREP_TIMING
  tim.change("cleanup");
#endif

  /* Begin loop over generalized contractions. */
  ogc1 = 0;
  for (i=0; i<nc1; i++) {
    tam1 = int_shell1->am(i) + dam1;
    if (tam1 < 0) continue;
    int tsize1 = INT_NCART_NN(tam1);
    ogc2 = 0;
    for (j=0; j<nc2; j++) {
      tam2 = int_shell2->am(j) + dam2;
      if (tam2 < 0) continue;
      int tsize2 = INT_NCART_NN(tam2);
      ogc3 = 0;
      for (k=0; k<nc3; k++) {
        tam3 = int_shell3->am(k) + dam3;
        if (tam3 < 0) continue;
        int tsize3 = INT_NCART_NN(tam3);
        ogc4 = 0;
        for (l=0; l<nc4; l++) {
          tam4 = int_shell4->am(l) + dam4;
          if (tam4 < 0) continue;
          int tsize4 = INT_NCART_NN(tam4);

#ifdef EREP_TIMING
  tim.change("shift");
#endif
  /* Shift angular momentum from 1 to 2 and from 3 to 4. */
  double *shiftbuffer = int_shiftgcam(i,j,k,l,tam1,tam2,tam3,tam4, eAB);
#ifdef EREP_TIMING
  tim.change("cleanup");
#endif

  /* Place the integrals in the integral buffer. */
  /* If permute_ is not set, then repack the integrals while copying. */
  if ((!permute_)&&(p12||p34||p13p24)) {
    int pam1,pam2,pam3,pam4;
    int psize234,psize34;
    int pogc1,pogc2,pogc3,pogc4;
    int psize1,psize2,psize3,psize4;
    pam1 = tam1;
    pam2 = tam2;
    pam3 = tam3;
    pam4 = tam4;
    pogc1 = ogc1;
    pogc2 = ogc2;
    pogc3 = ogc3;
    pogc4 = ogc4;
    psize1 = size1;
    psize2 = size2;
    psize3 = size3;
    psize4 = size4;
    if (p13p24) {
      iswtch(&pam1,&pam3);
      iswtch(&pam2,&pam4);
      iswtch(&pogc1,&pogc3);
      iswtch(&pogc2,&pogc4);
      iswtch(&psize1,&psize3);
      iswtch(&psize2,&psize4);
      }
    if (p34) {
      iswtch(&pam3,&pam4);
      iswtch(&pogc3,&pogc4);
      iswtch(&psize3,&psize4);
      }
    if (p12) {
      iswtch(&pam1,&pam2);
      iswtch(&pogc1,&pogc2);
      iswtch(&psize1,&psize2);
      }
    psize34 = psize4 * psize3;
    psize234 = psize34 * psize2;
    redundant_index = 0;
    int newindexoffset = pogc1*psize234 + pogc2*psize34 + pogc3*psize4 + pogc4;
    if (p13p24||p34) {
      int stride1=psize234;
      int stride2=psize34;
      int stride3=psize4;
      int stride4=1;
      int tmp;
      if (p12) {
        tmp=stride1; stride1=stride2; stride2=tmp;
        }
      if (p34) {
        tmp=stride3; stride3=stride4; stride4=tmp;
        }
      if (p13p24) {
        tmp=stride1; stride1=stride3; stride3=tmp;
        tmp=stride2; stride2=stride4; stride4=tmp;
        }
      int newindex1 = newindexoffset;
      for (int ci1=0; ci1<tsize1; ci1++) {
        int newindex12 = newindex1;
        for (int ci2=0; ci2<tsize2; ci2++) {
          int newindex123 = newindex12;
          for (int ci3=0; ci3<tsize3; ci3++) {
            double *tmp_shiftbuffer = &shiftbuffer[redundant_index];
            int newindex1234 = newindex123;
            for (int ci4=0; ci4<tsize4; ci4++) {
              int_buffer[newindex1234] = tmp_shiftbuffer[ci4];
              newindex1234 += stride4;
              }
            redundant_index+=tsize4;
            newindex123 += stride3;
            }
          newindex12 += stride2;
          }
        newindex1 += stride1;
        }
      }
    else if (nc3 == 1 && nc4 == 1) {
      // this is the p12 only case w/o gen contractions on 3 & 4
      // this special case collapses the 3rd and 4th indices together
      for (int ci1=0; ci1<tsize1; ci1++) {
        for (int ci2=0; ci2<tsize2; ci2++) {
          int pci1=ci1;
          int pci2=ci2;
          if (p12) {
            int tmp=pci1; pci1=pci2; pci2=tmp;
            }
          int newindex123 = newindexoffset + pci1*psize234 + pci2*psize34;
          double *tmp_int_buffer = &int_buffer[newindex123];
          double *tmp_shiftbuffer = &shiftbuffer[redundant_index];
          for (int ci34=0; ci34<psize34; ci34++)
            tmp_int_buffer[ci34] = tmp_shiftbuffer[ci34];
          redundant_index += psize34;
          }
        }
      }
    else {
      // this is the p12 only case w/gen. contr. on 3 & 4
      for (int ci1=0; ci1<tsize1; ci1++) {
        for (int ci2=0; ci2<tsize2; ci2++) {
          int pci1=ci1;
          int pci2=ci2;
          if (p12) {
            int tmp=pci1; pci1=pci2; pci2=tmp;
            }
          int newindex123 = newindexoffset + pci1*psize234 + pci2*psize34;
          for (int ci3=0; ci3<tsize3; ci3++) {
            double *tmp_int_buffer = &int_buffer[newindex123];
            double *tmp_shiftbuffer = &shiftbuffer[redundant_index];
            for (int ci4=0; ci4<tsize4; ci4++) {
              tmp_int_buffer[ci4] = tmp_shiftbuffer[ci4];
              }
            redundant_index += tsize4;
            newindex123 += psize4;
            }
          }
        }
      }
    }
  else if (nc3 == 1 && nc4 == 1) {
    // this special case collapses the 3rd and 4th indices together
    int size34 =  size3 * size4;
    int size234 = size2 * size34;
    double* redund_ints = shiftbuffer;
    redundant_index = 0;
    int newindex1 = ogc1*size234 + ogc2*size34 + ogc3*size4 + ogc4;
    for (int ci1=0; ci1<tsize1; ci1++) {
      int newindex12 = newindex1;
      for (int ci2=0; ci2<tsize2; ci2++) {
        double *tmp_int_buffer = &int_buffer[newindex12];
        double *tmp_redund_ints = &redund_ints[redundant_index];
        for (int ci34=0; ci34<size34; ci34++)
          tmp_int_buffer[ci34] = tmp_redund_ints[ci34];
        redundant_index += size34;
        newindex12 += size34;
        }
      newindex1 += size234;
      }
    }
  else {
    int size34 =  size3 * size4;
    int size234 = size2 * size34;
    double* redund_ints = shiftbuffer;
    redundant_index = 0;
    int newindex1 = ogc1*size234 + ogc2*size34 + ogc3*size4 + ogc4;
    for (int ci1=0; ci1<tsize1; ci1++) {
      int newindex12 = newindex1;
      for (int ci2=0; ci2<tsize2; ci2++) {
        int newindex123 = newindex12;
        for (int ci3=0; ci3<tsize3; ci3++) {
          double *tmp_int_buffer = &int_buffer[newindex123];
          double *tmp_redund_ints = &redund_ints[redundant_index];
          for (int ci4=0; ci4<tsize4; ci4++) {
            tmp_int_buffer[ci4] = tmp_redund_ints[ci4];
            }
          redundant_index += tsize4;
          newindex123 += size4;
          }
        newindex12 += size34;
        }
      newindex1 += size234;
      }
    }

    /* End loop over generalized contractions. */
          ogc4 += tsize4;
          }
        ogc3 += tsize3;
        }
      ogc2 += tsize2;
      }
    ogc1 += tsize1;
    }

  if (   !int_unit2
      && !int_unit4
      && int_integral_storage) {
#ifdef EREP_TIMING
      tim.change("maybe store");
#endif
      int_store_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24,size);
    }

  /* We branch here if an integral was precomputed and the int_buffer
   * has been already filled. */
  post_computation:

#ifdef EREP_TIMING
  tim.change("post");
#endif

  /* Unpermute all of the permuted quantities. */
  if ((!permute_)&&(p12||p34||p13p24)) {
    if (p13p24) {
      iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
      iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
      iswtch(&int_unit2,&int_unit4);
      iswtch(&am1,&am3);
      iswtch(&am2,&am4);
      iswtch(&am12,&am34);
      sswtch(&int_shell1,&int_shell3);
      swtch(pbs1,pbs3);
      sswtch(&int_shell2,&int_shell4);
      swtch(pbs2,pbs4);
      iswtch(&int_expweight1,&int_expweight3);
      iswtch(&int_expweight2,&int_expweight4);
      }
    if (p34) {
      iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&am3,&am4);
      sswtch(&int_shell3,&int_shell4);
      swtch(pbs3,pbs4);
      iswtch(&int_expweight3,&int_expweight4);
      }
    if (p12) {
      iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
      iswtch(&am1,&am2);
      sswtch(&int_shell1,&int_shell2);
      swtch(pbs1,pbs2);
      iswtch(&int_expweight1,&int_expweight2);
      }
    }

  pbs1_ = 0;
  pbs2_ = 0;
  pbs3_ = 0;
  pbs4_ = 0;

  /* Transform to pure am (if requested in the centers structure). */
  if (!(flags&INT_NOPURE)) {
      transform_2e(integral_, int_buffer, int_buffer,
                   &bs1_->shell(sh1),
                   int_unit2?int_unit_shell:&bs2_->shell(sh2),
                   &bs3_->shell(sh3),
                   int_unit4?int_unit_shell:&bs4_->shell(sh4));
    }

  /* Remove the redundant integrals, unless redundant_ is set. */
  if (!redundant_) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    if ((osh1 == osh4)&&(osh2 == osh3)&&(osh1 != osh2)) {
      ExEnv::errn() << scprintf("nonredundant integrals cannot be generated\n");
      fail();
      }
    e12 = (int_unit2?0:(osh1 == osh2));
    e13e24 = ((osh1 == osh3)
              && ((int_unit2 && int_unit4)
                  || ((int_unit2||int_unit4)?0:(osh2 == osh4))));
    e34 = (int_unit4?0:(osh3 == osh4));
    if (e12||e34||e13e24) {
      nonredundant_erep(int_buffer,e12,e34,e13e24,
                        int_shell1->nfunction(),
                        int_shell2->nfunction(),
                        int_shell3->nfunction(),
                        int_shell4->nfunction(),
                        &redundant_offset,
                        &nonredundant_offset);
      }
    }
    
#ifdef EREP_TIMING
  tim.exit("post");
  tim.exit(section);
#endif
  }

/* This computes the two electron derivatives for all unique
 * centers in the passed shell quartet.  One center in
 * the set of unique centers is not included.  This can
 * be computed as minus the sum of the other derivatives.
 * The list of centers for which integrals were computed can
 * be determined from the contents of der_centers.
 * The results are put into the global integral buffer in the
 * format:
 * +------------------+
 * | dercenter1 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 * | dercenter2 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 * | dercenter3 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 */

void
Int2eV3::erep_all1der(int &psh1, int &psh2, int &psh3, int &psh4,
                      der_centersv3_t *der_centers)
{
  double *current_buffer;
  int nints;
  double *user_int_buffer;
  int omit;
  GaussianBasisSet *cs[4];
  int sh[4];
  int n_unique;
  int i,j;
  GaussianShell *shell1,*shell2,*shell3,*shell4;
  GaussianBasisSet *ucs[4]; /* The centers struct for the unique centers. */
  int unum[4];        /* The number of times that this unique center occurs. */
  int uam[4];         /* The total angular momentum on each unique center. */
  int am[4];
  int osh[4];
  int cen[4];
  int ucen[4];
  int ncart;
  double *current_pure_buffer;

  cs[0] = bs1_.pointer();
  cs[1] = bs2_.pointer();
  cs[2] = bs3_.pointer();
  cs[3] = bs4_.pointer();

  sh[0] = psh1;
  sh[1] = psh2;
  sh[2] = psh3;
  sh[3] = psh4;

  /* Set up pointers to the current shells. */
  shell1 = &bs1_->shell(psh1);
  shell2 = &bs2_->shell(psh2);
  shell3 = &bs3_->shell(psh3);
  shell4 = &bs4_->shell(psh4);

  /* Number of cartesian and pure integrals. */
  ncart = shell1->ncartesian()*shell2->ncartesian()
         *shell3->ncartesian()*shell4->ncartesian();
  nints = shell1->nfunction()*shell2->nfunction()
         *shell3->nfunction()*shell4->nfunction();

  am[0] = shell1->max_am();
  am[1] = shell2->max_am();
  am[2] = shell3->max_am();
  am[3] = shell4->max_am();

  /* Compute the offset shell numbers. */
  osh[0] = psh1 + bs1_shell_offset_;
  osh[1] = psh2 + bs2_shell_offset_;
  osh[2] = psh3 + bs3_shell_offset_;
  osh[3] = psh4 + bs4_shell_offset_;

  for (i=0; i<4; i++) cen[i] = cs[i]->shell_to_center(sh[i]);

  /* This macro returns true if two shell centers are the same. */
#define SC(cs1,shc1,cs2,shc2) (((cs1)==(cs2))&&((shc1)==(shc2)))

  /* Build the list of unique centers structures and shells. */
  n_unique = 0;
  for (i=0; i<4; i++) {
    int unique = 1;
    for (j=0; j<n_unique; j++) {
      if (SC(ucs[j],ucen[j],cs[i],cen[i])) {
        unique = 0;
        uam[j] += am[i];
        unum[j]++;
        break;
        }
      }
    if (unique) {
      ucs[n_unique] = cs[i];
      ucen[n_unique] = cen[i];
      uam[n_unique] = am[i];
      unum[n_unique] = 1;
      n_unique++;
      }
    }

  /* Choose which unique center is to be omitted from the calculation. */
  if (n_unique == 1) {
    omit = 0;
    }
  else if (n_unique == 2) {
    if (unum[0]==3) omit = 0;
    else if (unum[1]==3) omit = 1;
    else if (uam[1]>uam[0]) omit = 1;
    else omit = 0;
    }
  else if (n_unique == 3) {
    if (unum[0]==2) omit = 0;
    else if (unum[1]==2) omit = 1;
    else omit = 2;
    }
  else {
    int max = 0;
    omit = 0;
    for (i=0; i<4; i++) {
      if (uam[i]>max) {
        max = uam[i];
        omit = i;
        }
      }
    }

  /* Save the location of the int_buffer. */
  user_int_buffer = int_buffer;
  int_buffer = int_derint_buffer;

  /* Zero out the result integrals. */
  for (i=0; i<3*(n_unique-1)*ncart; i++) user_int_buffer[i] = 0.0;

  /* Loop thru the unique centers, computing the integrals and
   * skip the derivative on the unique center specified by omit. */
  der_centers->n = 0;
  current_buffer = user_int_buffer;
  for (i=0; i<n_unique; i++) {
    if (i==omit) continue;

    der_centers->cs[der_centers->n] = ucs[i];
    der_centers->num[der_centers->n] = ucen[i];
    der_centers->n++;

    for (j=0; j<4; j++) {
      if (SC(ucs[i],ucen[i],cs[j],cen[j])) {
        int old_perm = permute();
        set_permute(0);
        compute_erep_1der(0,current_buffer,
                          &psh1,&psh2,&psh3,&psh4,j);
        set_permute(old_perm);
        }
      }

    current_buffer = &current_buffer[3*ncart];
    }

  /* Put the information about the omitted center into der_centers. */
  der_centers->ocs = ucs[omit];
  der_centers->onum = ucen[omit];

  /* Transform to pure am. */
  current_buffer = user_int_buffer;
  current_pure_buffer = user_int_buffer;
  for (i=0; i<3*der_centers->n; i++) {
      transform_2e(integral_, current_buffer, current_pure_buffer,
                   shell1, shell2, shell3, shell4);
      current_buffer = &current_buffer[ncart];
      current_pure_buffer = &current_pure_buffer[nints];
    }

  /* Eliminate redundant integrals, unless flags specifies otherwise. */
  current_buffer = user_int_buffer;
  if (!redundant_) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    int e12,e13e24,e34;
    int i;

    if ((osh[0] == osh[3])&&(osh[1] == osh[2])&&(osh[0] != osh[1])) {
      ExEnv::errn() << scprintf("nonredundant integrals cannot be generated (1der)\n");
      fail();
      }

    /* Shell equivalence information. */
    e12 = (osh[0] == osh[1]);
    e13e24 = ((osh[0] == osh[2]) && (osh[1] == osh[3]));
    e34 = (osh[2] == osh[3]);
    if (e12||e13e24||e34) {
      /* Repack x, y, and z integrals. */
      for (i=0; i<3*der_centers->n; i++) {
        nonredundant_erep(current_buffer,e12,e34,e13e24,
                          shell1->nfunction(),
                          shell2->nfunction(),
                          shell3->nfunction(),
                          shell4->nfunction(),
                          &redundant_offset,
                          &nonredundant_offset);
        }
      }
    }

  /* Return the integral buffers to their normal state. */
  int_derint_buffer = int_buffer;
  int_buffer = user_int_buffer;
  }

/* This will call compute_erep
 * to compute the derivatives in terms of order == 0 integrals.
 * flags are the flags to the int_comp_erep routine
 * psh1-4 are pointers to the shell numbers
 * dercenter is 0, 1, 2, or 3 -- the center within the integral
 *           which we are taking the derivative with respect to.
 * The results are accumulated in buffer, which cannot be the same
 * as the current int_buffer.
 */
void
Int2eV3::compute_erep_1der(int flags, double *buffer,
                           int *psh1, int *psh2, int *psh3, int *psh4,
                           int dercenter)
{
  int ii;
  int index;
  int size1,size2,size3,size4,size1234;
  int sizem234,sizem34,sizem2,sizem3,sizem4;
  int sizep234,sizep34,sizep2,sizep3,sizep4;
  GaussianShell *shell1,*shell2,*shell3,*shell4;

#ifdef DER_TIMING
  Timer tim("erep_1der");
#endif

  /* Set up pointers to the current shells. */
  shell1 = &bs1_->shell(*psh1);
  shell2 = &bs2_->shell(*psh2);
  shell3 = &bs3_->shell(*psh3);
  shell4 = &bs4_->shell(*psh4);

  if ((dercenter<0) || (dercenter > 3)) {
    ExEnv::errn() << scprintf("illegal derivative center -- must be 0, 1, 2, or 3\n");
    fail();
    }

  /* Offsets for the intermediates with original angular momentum. */
  for (ii=size1=0; ii<shell1->ncontraction(); ii++)
    size1 += INT_NCART(shell1->am(ii));
  for (ii=size2=0; ii<shell2->ncontraction(); ii++)
    size2 += INT_NCART(shell2->am(ii));
  for (ii=size3=0; ii<shell3->ncontraction(); ii++)
    size3 += INT_NCART(shell3->am(ii));
  for (ii=size4=0; ii<shell4->ncontraction(); ii++)
    size4 += INT_NCART(shell4->am(ii));

  size1234 = size1*size2*size3*size4;

#define DCTEST(n) ((dercenter==n)?1:0)
  /* Offsets for the intermediates with angular momentum decremented. */
  for (ii=sizem2=0; ii<shell2->ncontraction(); ii++) 
    sizem2 += INT_NCART(shell2->am(ii)-DCTEST(1));
  for (ii=sizem3=0; ii<shell3->ncontraction(); ii++) 
    sizem3 += INT_NCART(shell3->am(ii)-DCTEST(2));
  for (ii=sizem4=0; ii<shell4->ncontraction(); ii++) 
    sizem4 += INT_NCART(shell4->am(ii)-DCTEST(3));
  sizem34 = sizem4 * sizem3;
  sizem234 = sizem34 * sizem2;

  /* Offsets for the intermediates with angular momentum incremented. */
  for (ii=sizep2=0; ii<shell2->ncontraction(); ii++) 
    sizep2 += INT_NCART(shell2->am(ii)+DCTEST(1));
  for (ii=sizep3=0; ii<shell3->ncontraction(); ii++) 
    sizep3 += INT_NCART(shell3->am(ii)+DCTEST(2));
  for (ii=sizep4=0; ii<shell4->ncontraction(); ii++) 
    sizep4 += INT_NCART(shell4->am(ii)+DCTEST(3));
  sizep34 = sizep4 * sizep3;
  sizep234 = sizep34 * sizep2;

#ifdef DER_TIMING
  tim.enter("- erep");
#endif

  int old_perm = permute();
  set_permute(0);
  int old_red = redundant();
  set_redundant(1);
  compute_erep(flags|INT_NOPURE,
               psh1,psh2,psh3,psh4,
                   -DCTEST(0),
                   -DCTEST(1),
                   -DCTEST(2),
                   -DCTEST(3));
  set_permute(old_perm);
  set_redundant(old_red);

  /* Trouble if cpp is nonANSI. */
#define DERLOOP(index,indexp1) \
   oc##indexp1 = 0;\
   for ( c##indexp1 =0; c##indexp1 <shell##indexp1->ncontraction(); c##indexp1 ++) {\
     am[index] = shell##indexp1->am(c##indexp1);\
     FOR_CART(i[index],j[index],k[index],am[index])

#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCTEST(index));\
     }

#define ALLDERLOOPS\
    DERLOOP(0,1)\
      DERLOOP(1,2)\
        DERLOOP(2,3)\
          DERLOOP(3,4)

#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

  /* Place the contributions into the user integral buffer. */
  index = 0;

  if (dercenter==0) {
    int ogc1,ogc1m,gc1,i1,k1,f234,size234;
    size234=size2*size3*size4;

#ifdef DER_TIMING
  tim.change("- 0");
#endif
    /* The center 0 d/dx integrals */
    ogc1 = 0;
    ogc1m = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      // only integrals with x^n n>0 on center 0 contribute
      // so skip over n==0
      index += (am1+1)*size234;
      int c1 = am1+1;
      for (i1=1; i1<=am1; i1++) {
        double factor = -i1;
        for (k1=0; k1<=am1-i1; k1++) {
          int c1xm1 = c1-am1-1;//=INT_CARTINDEX(am1-1,i1-1,j1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1m+c1xm1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += factor * tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        }
      ogc1 += c1;
      ogc1m += INT_NCART(am1-1);
      }

    /* The center 0 d/dy integrals */
    ogc1 = 0;
    ogc1m = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      // only integrals with y^n n>0 on center 0 contribute
      // so skip over n==0 (by making i1+k1<am1)
      int c1 = 0;
      for (i1=0; i1<=am1; i1++) {
        for (k1=0; k1<=am1-i1-1; k1++) {
          double factor = -(am1-i1-k1);
          int c1ym1 = c1-i1;//=INT_CARTINDEX(am1-1,i1,j1-1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1m+c1ym1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += factor * tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        // account for the y^n n==0 case by an extra increment
        c1++;
        index+=size234;
        }
      ogc1 += c1;
      ogc1m += INT_NCART(am1-1);
      }

    /* The center 0 d/dz integrals */
    ogc1 = 0;
    ogc1m = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      int c1 = 0;
      for (i1=0; i1<=am1; i1++) {
        // only integrals with z^n n>0 on center 0 contribute
        // so skip over n==0
        c1++;
        index+=size234;
        for (k1=1; k1<=am1-i1; k1++) {
          double factor = -k1;
          int c1zm1 = c1-i1-1;//=INT_CARTINDEX(am1-1,i1,j1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1m+c1zm1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += factor * tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        }
      ogc1 += c1;
      ogc1m += INT_NCART(am1-1);
      }
    }
  else if (dercenter == 1) {
    int ogc2,ogc2m,gc2,i2,k2,f1,f34,size34,size234;
    size34 = size3*size4;
    size234 = size2*size3*size4;

#ifdef DER_TIMING
  tim.change("- 1");
#endif
    /* The center 1 d/dx integrals */
    ogc2 = 0;
    ogc2m = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      // only integrals with x^n n>0 on center 1 contribute
      // so skip over n==0
      int c2 = am2+1;
      for (i2=1; i2<=am2; i2++) {
        double factor = -i2;
        for (k2=0; k2<=am2-i2; k2++) {
          int c2xm1 = c2-am2-1;//=INT_CARTINDEX(am2-1,i2-1,j2)
          int buffer_index = (ogc2+c2)*size34;
          int int_buffer_index = (ogc2m+c2xm1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += factor * tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizem234;
            }
          c2++;
          }
        }
      ogc2 += c2;
      ogc2m += INT_NCART(am2-1);
      }
    index += size1234;

     /* The center 1 d/dy integrals */
    ogc2 = 0;
    ogc2m = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      // only integrals with y^n n>0 on center 1 contribute
      // so skip over n==0
      int c2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2-1; k2++) {
          double factor = -(am2-k2-i2);
          int c2ym1 = c2-i2;//=INT_CARTINDEX(am2-1,i2,j2-1)
          int buffer_index = size1234 + (ogc2+c2)*size34;
          int int_buffer_index = (ogc2m+c2ym1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += factor * tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizem234;
            }
          c2++;
          }
        // account for the y^n n==0 case by an extra increment
        c2++;
        }
      ogc2 += c2;
      ogc2m += INT_NCART(am2-1);
      }
    index += size1234;

     /* The center 1 d/dz integrals */
    ogc2 = 0;
    ogc2m = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      // only integrals with z^n n>0 on center 1 contribute
      // so skip over n==0
      int c2 = 0;
      for (i2=0; i2<=am2; i2++) {
        // account for the z^n n==0 case by an extra increment
        c2++;
        for (k2=1; k2<=am2-i2; k2++) {
          double factor = -k2;
          int c2zm1 = c2-i2-1;//=INT_CARTINDEX(am2-1,i2,j2-1)
          int buffer_index = size1234+size1234+(ogc2+c2)*size34;
          int int_buffer_index = (ogc2m+c2zm1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += factor * tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizem234;
            }
          c2++;
          }
        }
      ogc2 += c2;
      ogc2m += INT_NCART(am2-1);
      }
    index += size1234;
    }
  else if (dercenter == 2) {
    int ogc3,ogc3m,gc3,i3,k3,f12,f4,size12,size34;
    size12 = size1*size2;
    size34 = size3*size4;

#ifdef DER_TIMING
  tim.change("- 2");
#endif
    /* The center 2 d/dx integrals */
    ogc3 = 0;
    ogc3m = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      // only integrals with x^n n>0 on center 2 contribute
      // so skip over n==0
      int c3 = am3+1;
      for (i3=1; i3<=am3; i3++) {
        double factor = -i3;
        for (k3=0; k3<=am3-i3; k3++) {
          int c3xm1 = c3-am3-1;//=INT_CARTINDEX(am3-1,i3-1,j3)
          int buffer_index = (ogc3+c3)*size4;
          int int_buffer_index = (ogc3m+c3xm1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += factor * tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizem34;
            }
          c3++;
          }
        }
      ogc3 += c3;
      ogc3m += INT_NCART(am3-1);
      }
    index += size1234;

     /* The center 2 d/dy integrals */
    ogc3 = 0;
    ogc3m = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      // only integrals with y^n n>0 on center 2 contribute
      // so skip over n==0
      int c3 = 0;
      for (i3=0; i3<=am3; i3++) {
        for (k3=0; k3<=am3-i3-1; k3++) {
          double factor = -(am3-k3-i3);
          int c3ym1 = c3-i3;//=INT_CARTINDEX(am3-1,i3,j3-1)
          int buffer_index = size1234 + (ogc3+c3)*size4;
          int int_buffer_index = (ogc3m+c3ym1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += factor * tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizem34;
            }
          c3++;
          }
        // account for the y^n n==0 case by an extra increment
        c3++;
        }
      ogc3 += c3;
      ogc3m += INT_NCART(am3-1);
      }
    index += size1234;

     /* The center 2 d/dz integrals */
    ogc3 = 0;
    ogc3m = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      // only integrals with z^n n>0 on center 2 contribute
      // so skip over n==0
      int c3 = 0;
      for (i3=0; i3<=am3; i3++) {
        // account for the z^n n==0 case by an extra increment
        c3++;
        for (k3=1; k3<=am3-i3; k3++) {
          double factor = -k3;
          int c3zm1 = c3-i3-1;//=INT_CARTINDEX(am3-1,i3,j3)
          int buffer_index = size1234+size1234+(ogc3+c3)*size4;
          int int_buffer_index = (ogc3m+c3zm1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += factor * tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizem34;
            }
          c3++;
          }
        }
      ogc3 += c3;
      ogc3m += INT_NCART(am3-1);
      }
    index += size1234;
    }
  else if (dercenter == 3) {
    int ogc4,ogc4m,gc4,i4,k4,f123,size123;
    size123 = size1*size2*size3;

#ifdef DER_TIMING
  tim.change("- 3");
#endif
    /* The center 3 d/dx integrals */
    ogc4 = 0;
    ogc4m = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      // only integrals with x^n n>0 on center 3 contribute
      // so skip over n==0
      int c4 = am4+1;
      for (i4=1; i4<=am4; i4++) {
        double factor = -i4;
        for (k4=0; k4<=am4-i4; k4++) {
          int c4xm1 = c4-am4-1;//=INT_CARTINDEX(am4-1,i4-1,j4)
          int buffer_index = ogc4+c4;
          int int_buffer_index = ogc4m+c4xm1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += factor * int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizem4;
            }
          c4++;
          }
        }
      ogc4 += c4;
      ogc4m += INT_NCART(am4-1);
      }
    index += size1234;

     /* The center 3 d/dy integrals */
    ogc4 = 0;
    ogc4m = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      // only integrals with y^n n>0 on center 3 contribute
      // so skip over n==0
      int c4 = 0;
      for (i4=0; i4<=am4; i4++) {
        for (k4=0; k4<=am4-i4-1; k4++) {
          double factor = -(am4-k4-i4);
          int c4ym1 = c4-i4;//=INT_CARTINDEX(am4-1,i4,j4-1)
          int buffer_index = size1234 + ogc4+c4;
          int int_buffer_index = ogc4m+c4ym1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += factor * int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizem4;
            }
          c4++;
          }
        // account for the y^n n==0 case by an extra increment
        c4++;
        }
      ogc4 += c4;
      ogc4m += INT_NCART(am4-1);
      }
    index += size1234;

    /* The center 3 d/dz integrals */
    ogc4 = 0;
    ogc4m = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      // only integrals with z^n n>0 on center 3 contribute
      // so skip over n==0
      int c4 = 0;
      for (i4=0; i4<=am4; i4++) {
        // account for the z^n n==0 case by an extra increment
        c4++;
        for (k4=1; k4<=am4-i4; k4++) {
          double factor = -k4;
          int c4zm1 = c4-i4-1;//=INT_CARTINDEX(am4-1,i4,j4-1)
          int buffer_index = size1234+size1234+ogc4+c4;
          int int_buffer_index = ogc4m+c4zm1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += factor * int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizem4;
            }
          c4++;
          }
        }
      ogc4 += c4;
      ogc4m += INT_NCART(am4-1);
      }
    index += size1234;
    }

#ifdef DER_TIMING
  tim.change("+ erep");
#endif

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  if (dercenter==0) int_expweight1 = 1;
  else if (dercenter==1) int_expweight2 = 1;
  else if (dercenter==2) int_expweight3 = 1;
  else if (dercenter==3) int_expweight4 = 1;
  old_perm = permute();
  set_permute(0);
  old_red = redundant();
  set_redundant(1);
  compute_erep(flags|INT_NOPURE,
               psh1,psh2,psh3,psh4,
                     DCTEST(0),
                     DCTEST(1),
                     DCTEST(2),
                     DCTEST(3));
  set_permute(old_perm);
  set_redundant(old_red);
  if (dercenter==0) int_expweight1 = 0;
  else if (dercenter==1) int_expweight2 = 0;
  else if (dercenter==2) int_expweight3 = 0;
  else if (dercenter==3) int_expweight4 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  if (dercenter==0) {
    int ogc1,ogc1p,gc1,i1,k1,f234,size234;
    size234=size2*size3*size4;

#ifdef DER_TIMING
  tim.change("+ 0");
#endif
    /* The center 0 d/dx integrals */
    ogc1 = 0;
    ogc1p = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      int c1 = 0;
      for (i1=0; i1<=am1; i1++) {
        for (k1=0; k1<=am1-i1; k1++) {
          int c1xp1 = c1+am1+2;//=INT_CARTINDEX(am1+1,i1+1,j1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1p+c1xp1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        }
      ogc1 += c1;
      ogc1p += INT_NCART(am1+1);
      }

    /* The center 0 d/dy integrals */
    ogc1 = 0;
    ogc1p = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      int c1 = 0;
      for (i1=0; i1<=am1; i1++) {
        for (k1=0; k1<=am1-i1; k1++) {
          int c1yp1 = c1+i1;//=INT_CARTINDEX(am1+1,i1,j1+1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1p+c1yp1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        }
      ogc1 += c1;
      ogc1p += INT_NCART(am1+1);
      }

    /* The center 0 d/dz integrals */
    ogc1 = 0;
    ogc1p = 0;
    for (gc1=0; gc1<shell1->ncontraction(); gc1++) {
      int am1 = shell1->am(gc1);
      int c1 = 0;
      for (i1=0; i1<=am1; i1++) {
        for (k1=0; k1<=am1-i1; k1++) {
          int c1zp1 = c1+i1+1;//=INT_CARTINDEX(am1+1,i1,j1)
          double *tmp_buffer=&buffer[index];
          double *tmp_int_buffer=&int_buffer[(ogc1p+c1zp1)*size234];
          for (f234=0; f234<size234; f234++) {
            tmp_buffer[f234] += tmp_int_buffer[f234];
            }
          index+=size234;
          c1++;
          }
        }
      ogc1 += c1;
      ogc1p += INT_NCART(am1+1);
      }
    }
  else if (dercenter == 1) {
    int ogc2,ogc2p,gc2,i2,k2,f1,f34,size34,size234;
    size34 = size3*size4;
    size234 = size2*size3*size4;

#ifdef DER_TIMING
  tim.change("+ 1");
#endif
    /* The center 1 d/dx integrals */
    ogc2 = 0;
    ogc2p = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      int c2=0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          int c2xp1 = c2+am2+2;//=INT_CARTINDEX(am2+1,i2+1,j2)
          int buffer_index = (ogc2+c2)*size34;
          int int_buffer_index = (ogc2p+c2xp1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizep234;
            }
          c2++;
          }
        }
      ogc2 += c2;
      ogc2p += INT_NCART(am2+1);
      }
    index += size1234;

     /* The center 1 d/dy integrals */
    ogc2 = 0;
    ogc2p = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      int c2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          int c2yp1 = c2+i2;//=INT_CARTINDEX(am2+1,i2,j2+1)
          int buffer_index = size1234 + (ogc2+c2)*size34;
          int int_buffer_index = (ogc2p+c2yp1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizep234;
            }
          c2++;
          }
        }
      ogc2 += c2;
      ogc2p += INT_NCART(am2+1);
      }
    index += size1234;

    /* The center 1 d/dz integrals */
    ogc2 = 0;
    ogc2p = 0;
    for (gc2=0; gc2<shell2->ncontraction(); gc2++) {
      int am2 = shell2->am(gc2);
      int c2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          int c2zp1 = c2+i2+1;//=INT_CARTINDEX(am2+1,i2,j2+1)
          int buffer_index = size1234+size1234+(ogc2+c2)*size34;
          int int_buffer_index = (ogc2p+c2zp1)*size34;
          for (f1=0; f1<size1; f1++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f34=0; f34<size34; f34++) {
                tmp_buffer[f34] += tmp_int_buffer[f34];
              }
            buffer_index += size234;
            int_buffer_index += sizep234;
            }
          c2++;
          }
        }
      ogc2 += c2;
      ogc2p += INT_NCART(am2+1);
      }
    index += size1234;
    }
  else if (dercenter == 2) {
    int ogc3,ogc3p,gc3,i3,k3,f12,f4,size12,size34;
    size12 = size1*size2;
    size34 = size3*size4;

#ifdef DER_TIMING
  tim.change("+ 2");
#endif
    /* The center 2 d/dx integrals */
    ogc3 = 0;
    ogc3p = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      int c3 = 0;
      for (i3=0; i3<=am3; i3++) {
        for (k3=0; k3<=am3-i3; k3++) {
          int c3xp1 = c3+am3+2;//=INT_CARTINDEX(am3+1,i3+1,j3)
          int buffer_index = (ogc3+c3)*size4;
          int int_buffer_index = (ogc3p+c3xp1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizep34;
            }
          c3++;
          }
        }
      ogc3 += c3;
      ogc3p += INT_NCART(am3+1);
      }
    index += size1234;

     /* The center 2 d/dy integrals */
    ogc3 = 0;
    ogc3p = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      int c3 = 0;
      for (i3=0; i3<=am3; i3++) {
        for (k3=0; k3<=am3-i3; k3++) {
          int c3yp1 = c3+i3;//=INT_CARTINDEX(am3+1,i3,j3+1)
          int buffer_index = size1234 + (ogc3+c3)*size4;
          int int_buffer_index = (ogc3p+c3yp1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizep34;
            }
          c3++;
          }
        }
      ogc3 += c3;
      ogc3p += INT_NCART(am3+1);
      }
    index += size1234;

     /* The center 2 d/dz integrals */
    ogc3 = 0;
    ogc3p = 0;
    for (gc3=0; gc3<shell3->ncontraction(); gc3++) {
      int am3 = shell3->am(gc3);
      int c3 = 0;
      for (i3=0; i3<=am3; i3++) {
        for (k3=0; k3<=am3-i3; k3++) {
          int c3zp1 = c3+i3+1;//=INT_CARTINDEX(am3+1,i3,j3)
          int buffer_index = size1234+size1234+(ogc3+c3)*size4;
          int int_buffer_index = (ogc3p+c3zp1)*size4;
          for (f12=0; f12<size12; f12++) {
            double *tmp_buffer=&buffer[buffer_index];
            double *tmp_int_buffer=&int_buffer[int_buffer_index];
            for (f4=0; f4<size4; f4++) {
                tmp_buffer[f4] += tmp_int_buffer[f4];
              }
            buffer_index += size34;
            int_buffer_index += sizep34;
            }
          c3++;
          }
        }
      ogc3 += c3;
      ogc3p += INT_NCART(am3+1);
      }
    index += size1234;
    }
  else if (dercenter == 3) {
    int ogc4,ogc4p,gc4,i4,k4,f123,size123;
    size123 = size1*size2*size3;

#ifdef DER_TIMING
  tim.change("+ 3");
#endif
    /* The center 3 d/dx integrals */
    ogc4 = 0;
    ogc4p = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      int c4 = 0;
      for (i4=0; i4<=am4; i4++) {
        for (k4=0; k4<=am4-i4; k4++) {
          int c4xp1 = c4+am4+2;//=INT_CARTINDEX(am4+1,i4+1,j4)
          int buffer_index = ogc4+c4;
          int int_buffer_index = ogc4p+c4xp1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizep4;
            }
          c4++;
          }
        }
      ogc4 += c4;
      ogc4p += INT_NCART(am4+1);
      }
    index += size1234;

     /* The center 3 d/dy integrals */
    ogc4 = 0;
    ogc4p = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      int c4 = 0;
      for (i4=0; i4<=am4; i4++) {
        for (k4=0; k4<=am4-i4; k4++) {
          int c4yp1 = c4+i4;//=INT_CARTINDEX(am4+1,i4,j4+1)
          int buffer_index = size1234 + ogc4+c4;
          int int_buffer_index = ogc4p+c4yp1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizep4;
            }
          c4++;
          }
        }
      ogc4 += c4;
      ogc4p += INT_NCART(am4+1);
      }
    index += size1234;

    /* The center 3 d/dz integrals */
    ogc4 = 0;
    ogc4p = 0;
    for (gc4=0; gc4<shell4->ncontraction(); gc4++) {
      int am4 = shell4->am(gc4);
      int c4 = 0;
      for (i4=0; i4<=am4; i4++) {
        for (k4=0; k4<=am4-i4; k4++) {
          int c4zp1 = c4+i4+1;//=INT_CARTINDEX(am4+1,i4,j4)
          int buffer_index = size1234+size1234+ogc4+c4;
          int int_buffer_index = ogc4p+c4zp1;
          for (f123=0; f123<size123; f123++) {
            buffer[buffer_index] += int_buffer[int_buffer_index];
            buffer_index += size4;
            int_buffer_index += sizep4;
            }
          c4++;
          }
        }
      ogc4 += c4;
      ogc4p += INT_NCART(am4+1);
      }
    index += size1234;
    }
#ifdef DER_TIMING
  tim.exit(0);
  tim.exit(0);
#endif
  }

void
Int2eV3::nonredundant_erep(double *buffer, int e12, int e34, int e13e24,
                           int n1, int n2, int n3, int n4,
                           int *red_off, int *nonred_off)
{
  int nonredundant_index;
  int i,j,k,l;

  double *redundant_ptr = &buffer[*red_off];
  double *nonredundant_ptr = &buffer[*nonred_off];

  nonredundant_index = 0;
  int n34 = n3*n4;
  for (i=0; i<n1; i++) {
    int jmax = INT_MAX2(e12,i,n2);
    for (j=0; j<=jmax; j++) {
      int kmax = INT_MAX3(e13e24,i,n3);
      for (k=0; k<=kmax; k++) {
        int lmax = INT_MAX4(e13e24,e34,i,j,k,n4);
        for (l=0; l<=lmax; l++) {
          nonredundant_ptr[l] = redundant_ptr[l];
          }
        redundant_ptr += n4;
        nonredundant_index += lmax+1;
        nonredundant_ptr += lmax+1;
        }
      redundant_ptr += (n3-(kmax+1))*n4;
      }
    redundant_ptr += (n2-(jmax+1))*n34;
    }
  *red_off += n1*n2*n34;
  *nonred_off += nonredundant_index;
  }

/* This is an alternate interface to int_erep_all1der.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL, and the dercenters pointer. */
void
Int2eV3::erep_all1der(int *shells, int  *sizes,
                      der_centersv3_t *dercenters)
{
  erep_all1der(shells[0],shells[1],shells[2],shells[3],
               dercenters);
  if (sizes) {
    sizes[0] = bs1_->shell(shells[0]).nfunction();
    sizes[1] = bs2_->shell(shells[1]).nfunction();
    sizes[2] = bs3_->shell(shells[2]).nfunction();
    sizes[3] = bs4_->shell(shells[3]).nfunction();
    }
  }

void
Int2eV3::int_erep_bound1der(int flags, int bsh1, int bsh2, int *size)
{
  double *current_buffer;
  int nints;
  double *user_int_buffer;
  int i;
  GaussianShell *shell1,*shell2,*shell3,*shell4;
  int osh[4];
  int sh1 = bsh1;
  int sh2 = bsh2;
  int sh3 = bsh1;
  int sh4 = bsh2;
  int *psh1 = &sh1;
  int *psh2 = &sh2;
  int *psh3 = &sh3;
  int *psh4 = &sh4;
  int ncart;
  double *current_pure_buffer;

  /* Set up pointers to the current shells. */
  shell1 = &bs1_->shell(*psh1);
  shell2 = &bs2_->shell(*psh2);
  shell3 = &bs3_->shell(*psh3);
  shell4 = &bs4_->shell(*psh4);

  /* Number of cartesian and pure integrals. */
  ncart = shell1->ncartesian()*shell2->ncartesian()
         *shell3->ncartesian()*shell4->ncartesian();
  nints = shell1->nfunction() * shell2->nfunction()
        * shell3->nfunction() * shell4->nfunction();

  /* Compute the offset shell numbers. */
  osh[0] = *psh1 + bs1_shell_offset_;
  osh[1] = *psh2 + bs2_shell_offset_;
  osh[2] = *psh3 + bs3_shell_offset_;
  osh[3] = *psh4 + bs4_shell_offset_;

  /* Save the location of the int_buffer. */
  user_int_buffer = int_buffer;
  int_buffer = int_derint_buffer;

  /* Zero out the result integrals. */
  for (i=0; i<3*ncart; i++) user_int_buffer[i] = 0.0;

  /* Set the size so it is available to the caller. */
  *size = nints * 3;

  current_buffer = user_int_buffer;
  int old_perm = permute();
  compute_erep_bound1der(flags|INT_NOPURE,current_buffer,
                          psh1,psh2,psh3,psh4);
  set_permute(old_perm);

  /* Transform to pure am. */
  current_buffer = user_int_buffer;
  current_pure_buffer = user_int_buffer;
  for (i=0; i<3; i++) {
      transform_2e(integral_, current_buffer, current_pure_buffer,
                   &bs1_->shell(sh1),
                   &bs2_->shell(sh2),
                   &bs3_->shell(sh3),
                   &bs4_->shell(sh4));
      current_buffer = &current_buffer[ncart];
      current_pure_buffer = &current_pure_buffer[nints];
    }

  /* Eliminate redundant integrals, unless flags specifies otherwise. */
  current_buffer = user_int_buffer;
  if (!redundant_) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    int e12,e13e24,e34;
    int i;

    if ((osh[0] == osh[3])&&(osh[1] == osh[2])&&(osh[0] != osh[1])) {
      ExEnv::errn() << scprintf("nonredundant integrals cannot be generated (1der)\n");
      fail();
      }

    /* Shell equivalence information. */
    e12 = (osh[0] == osh[1]);
    e13e24 = ((osh[0] == osh[2]) && (osh[1] == osh[3]));
    e34 = (osh[2] == osh[3]);
    /* Repack x, y, and z integrals. */
    for (i=0; i<3; i++) {
      nonredundant_erep(current_buffer,e12,e34,e13e24,
                             shell1->nfunction(),
                             shell2->nfunction(),
                             shell3->nfunction(),
                             shell4->nfunction(),
                             &redundant_offset,
                             &nonredundant_offset);
      }
    }

  /* Return the integral buffers to their normal state. */
  int_derint_buffer = int_buffer;
  int_buffer = user_int_buffer;
  }

/* This routine computes a quantity needed to compute the
 * derivative integral bounds.
 * It fills int_buffer with (sh1+i sh2|sh1+i sh2).
 */
void
Int2eV3::compute_erep_bound1der(int flags, double *buffer,
                                int *psh1, int *psh2, int *psh3, int *psh4)
{
  int oc1,oc2,oc3,oc4;
  int ii;
  int c1,c2,c3,c4;
  int i[4],j[4],k[4],am[4];
  int index;
  int sizem234,sizem34,sizem2,sizem3,sizem4;
  int sizep234,sizep34,sizep2,sizep3,sizep4;
  int sizepm234,sizepm34,sizepm2,sizepm3,sizepm4;
  GaussianShell *shell1,*shell2,*shell3,*shell4;

  /* Set up pointers to the current shells. */
  shell1 = &bs1_->shell(*psh1);
  shell2 = &bs2_->shell(*psh2);
  shell3 = &bs3_->shell(*psh3);
  shell4 = &bs4_->shell(*psh4);

#define DCT1(n) ((((n)==0)||((n)==2))?1:0)
#define DCT2(n) ((((n)==0)||((n)==2))?((((n)==0))?1:-1):0)
  /* Offsets for the intermediates with angular momentum decremented. */
  for (ii=sizem2=0; ii<shell2->ncontraction(); ii++) 
    sizem2 += INT_NCART(shell2->am(ii)-DCT1(1));
  for (ii=sizem3=0; ii<shell3->ncontraction(); ii++) 
    sizem3 += INT_NCART(shell3->am(ii)-DCT1(2));
  for (ii=sizem4=0; ii<shell4->ncontraction(); ii++) 
    sizem4 += INT_NCART(shell4->am(ii)-DCT1(3));
  sizem34 = sizem4 * sizem3;
  sizem234 = sizem34 * sizem2;

  /* Offsets for the intermediates with angular momentum incremented. */
  for (ii=sizep2=0; ii<shell2->ncontraction(); ii++) 
    sizep2 += INT_NCART(shell2->am(ii)+DCT1(1));
  for (ii=sizep3=0; ii<shell3->ncontraction(); ii++) 
    sizep3 += INT_NCART(shell3->am(ii)+DCT1(2));
  for (ii=sizep4=0; ii<shell4->ncontraction(); ii++) 
    sizep4 += INT_NCART(shell4->am(ii)+DCT1(3));
  sizep34 = sizep4 * sizep3;
  sizep234 = sizep34 * sizep2;

  /* Offsets for the intermediates with angular momentum incremented and
   * decremented. */
  for (ii=sizepm2=0; ii<shell2->ncontraction(); ii++) 
    sizepm2 += INT_NCART(shell2->am(ii)+DCT2(1));
  for (ii=sizepm3=0; ii<shell3->ncontraction(); ii++) 
    sizepm3 += INT_NCART(shell3->am(ii)+DCT2(2));
  for (ii=sizepm4=0; ii<shell4->ncontraction(); ii++) 
    sizepm4 += INT_NCART(shell4->am(ii)+DCT2(3));
  sizepm34 = sizepm4 * sizepm3;
  sizepm234 = sizepm34 * sizepm2;

  /* END_DERLOOP must be redefined here because it previously depended
   * on the DCTEST macro */
#undef END_DERLOOP
#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCT1(index));\
     }

#undef END_ALLDERLOOPS
#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

  int old_perm = permute();
  set_permute(0);
  int old_red = redundant();
  set_redundant(1);
  compute_erep(flags,psh1,psh2,psh3,psh4,
                   -DCT1(0),
                   -DCT1(1),
                   -DCT1(2),
                   -DCT1(3));
  set_permute(old_perm);
  set_redundant(old_red);

   /* Place the contributions into the user integral buffer. */
   index = 0;
   /* The d/dx integrals */
  ALLDERLOOPS
    if (i[0]>0 && i[2]>0) {
      buffer[index] += i[0] * i[2] * int_buffer[
        (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0]-DCT1(0),j[0])) * sizem234
       +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1]-DCT1(1),j[1])) * sizem34
       +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2]-DCT1(2),j[2])) * sizem4
       +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3]-DCT1(3),j[3]))
       ];
      }
    index++;
    END_ALLDERLOOPS(-)

   /* The d/dy integrals */
  ALLDERLOOPS
    if (j[0]>0 && j[2]>0) {
    buffer[index] += j[0] * j[2] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0],j[0]-DCT1(0))) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1],j[1]-DCT1(1))) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2],j[2]-DCT1(2))) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3],j[3]-DCT1(3)))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

   /* The d/dz integrals */
  ALLDERLOOPS
    if (k[0]>0 && k[2]>0) {
    buffer[index] += k[0] * k[2] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0],j[0])) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1],j[1])) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2],j[2])) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3],j[3]))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  int_expweight1 = 1;
  int_expweight3 = 1;
  old_perm = permute();
  set_permute(0);
  old_red = redundant();
  set_redundant(1);
  compute_erep(flags,psh1,psh2,psh3,psh4,
                     DCT1(0),
                     DCT1(1),
                     DCT1(2),
                     DCT1(3));
  set_permute(old_perm);
  set_redundant(old_red);
  int_expweight1 = 0;
  int_expweight3 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  /* The d/dx integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT1(0),i[0]+DCT1(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCT1(1),i[1]+DCT1(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCT1(2),i[2]+DCT1(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCT1(3),i[3]+DCT1(3),j[3]))
            ];
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT1(0),i[0],j[0]+DCT1(0)))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCT1(1),i[1],j[1]+DCT1(1)))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCT1(2),i[2],j[2]+DCT1(2)))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCT1(3),i[3],j[3]+DCT1(3)))
            ];
          index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCT1(0),i[0],j[0])) * sizep234
              +(oc2 + INT_CARTINDEX(am[1]+DCT1(1),i[1],j[1])) * sizep34
              +(oc3 + INT_CARTINDEX(am[2]+DCT1(2),i[2],j[2])) * sizep4
              +(oc4 + INT_CARTINDEX(am[3]+DCT1(3),i[3],j[3]))
              ];
          index++;
    END_ALLDERLOOPS(+)

  /* END_DERLOOP must be redefined here because it previously depended
   * on the DCT1 macro */
#undef END_DERLOOP
#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCT2(index));\
     }

#undef END_ALLDERLOOPS
#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  int_expweight1 = 1;
  old_perm = permute();
  set_permute(0);
  old_red = redundant();
  set_redundant(1);
  compute_erep(flags,psh1,psh2,psh3,psh4,
                     DCT2(0),
                     DCT2(1),
                     DCT2(2),
                     DCT2(3));
  set_permute(old_perm);
  set_redundant(old_red);
  int_expweight1 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  /* The d/dx integrals */
  ALLDERLOOPS
    if (i[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT2(0),i[0]+DCT2(0),j[0]))*sizepm234
            +(oc2+INT_CARTINDEX(am[1]+DCT2(1),i[1]+DCT2(1),j[1]))*sizepm34
            +(oc3+INT_CARTINDEX(am[2]+DCT2(2),i[2]+DCT2(2),j[2]))*sizepm4
            +(oc4+INT_CARTINDEX(am[3]+DCT2(3),i[3]+DCT2(3),j[3]))
            ];
       }
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
    if (j[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT2(0),i[0],j[0]+DCT2(0)))*sizepm234
            +(oc2+INT_CARTINDEX(am[1]+DCT2(1),i[1],j[1]+DCT2(1)))*sizepm34
            +(oc3+INT_CARTINDEX(am[2]+DCT2(2),i[2],j[2]+DCT2(2)))*sizepm4
            +(oc4+INT_CARTINDEX(am[3]+DCT2(3),i[3],j[3]+DCT2(3)))
            ];
       }
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
    if (k[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCT2(0),i[0],j[0])) * sizepm234
              +(oc2 + INT_CARTINDEX(am[1]+DCT2(1),i[1],j[1])) * sizepm34
              +(oc3 + INT_CARTINDEX(am[2]+DCT2(2),i[2],j[2])) * sizepm4
              +(oc4 + INT_CARTINDEX(am[3]+DCT2(3),i[3],j[3]))
              ];
       }
    index++;
    END_ALLDERLOOPS(+)
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
