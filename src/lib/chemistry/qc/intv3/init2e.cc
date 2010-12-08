//
// init2e.cc
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

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>

using namespace std;
using namespace sc;

static void
fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

/* Initialize the 2e integral computation routines.
 * storage = the amount of storage available in bytes
 * order = order of derivative, must be zero or one
 * cs1 = center structure for center 1
 * cs2 = center structure for center 2
 * cs3 = center structure for center 3
 * cs4 = center structure for center 4
 * The integrals which will be computed are (cs1 cs2|cs3 cs4).
 * This function returns the pointer to the buffer where the
 * integrals are stored.
 */
double *
Int2eV3::int_initialize_erep(size_t storage, int order,
                             const Ref<GaussianBasisSet> &cs1,
                             const Ref<GaussianBasisSet> &cs2,
                             const Ref<GaussianBasisSet> &cs3,
                             const Ref<GaussianBasisSet> &cs4)
{
  int nc1,nc2,nc3,nc4;
  int jmax,jmax1,jmax2,jmax3,jmax4;

  redundant_ = 1;
  permute_ = 0;

  int_unit_shell = 0;

  /* Reset the integral storage variables. */
  int_integral_storage = 0;
  used_storage_ = 0;

  /* Turn off exponent weighted contractions. */
  int_expweight1 = 0;
  int_expweight2 = 0;
  int_expweight3 = 0;
  int_expweight4 = 0;

  /* See if the order of derivative needed is allowed. */
  if (order > 1) {
    ExEnv::errn() << scprintf("int_initialize_erep cannot handle order>1, yet\n");
    }

  if (order > 0) {
    int_derivative_bounds = 1;
    }
  else {
    int_derivative_bounds = 0;
    }

  /* Put the center pointers into the global centers pointers. */
  int_cs1 = cs1;
  int_cs2 = cs2;
  int_cs3 = cs3;
  int_cs4 = cs4;

  /* Find the max angular momentum on each center. */
  jmax1 = cs1->max_angular_momentum();
  if (!int_unit2) jmax2 = cs2->max_angular_momentum();
  else jmax2 = 0;
  jmax3 = cs3->max_angular_momentum();
  if (!int_unit4) jmax4 = cs4->max_angular_momentum();
  else jmax4 = 0;

  /* Find the maximum number of contractions in a shell on each center. */
  nc1 = cs1->max_ncontraction();
  if (!int_unit2) nc2 = cs2->max_ncontraction();
  else nc2 = 1;
  nc3 = cs3->max_ncontraction();
  if (!int_unit4) nc4 = cs4->max_ncontraction();
  else nc4 = 1;

  /* Initialize the Fj(T) routine. */
  jmax = jmax1+jmax2+jmax3+jmax4;
  if (int_derivative_bounds) {
      fjt_ = new FJT(jmax + 2*order); /* The 2 is for bounds checking */
    }
  else {
      fjt_ = new FJT(jmax + order);
    }

  /* Initialize the build and shift routines. */
  int_init_buildgc(order,jmax1,jmax2,jmax3,jmax4,nc1,nc2,nc3,nc4);
  int_init_shiftgc(order,jmax1,jmax2,jmax3,jmax4);

  /* Allocate storage for the integral buffer. */
  int maxsize = cs1->max_ncartesian_in_shell()
                *(int_unit2?1:cs2->max_ncartesian_in_shell())
                *cs3->max_ncartesian_in_shell()
                *(int_unit4?1:cs4->max_ncartesian_in_shell());
  if (order==0) {
    int_buffer = allocate<double>(maxsize);
    int_derint_buffer = 0;
    }
  else if (order==1) {
    int nderint;
    nderint = cs1->max_ncartesian_in_shell(1)
             *(int_unit2?1:cs2->max_ncartesian_in_shell(1))
             *cs3->max_ncartesian_in_shell(1)
             *(int_unit4?1:cs4->max_ncartesian_in_shell(1));

    /* Allocate the integral buffers. */
    int_buffer = allocate<double>(9*maxsize);
    int_derint_buffer = allocate<double>(nderint);
    if (!int_derint_buffer) {
      ExEnv::errn() << scprintf("couldn't malloc intermed storage for derivative ints\n");
      fail();
      }
    }

  if (!int_buffer) {
    ExEnv::errn() << scprintf("couldn't allocate integrals\n");
    fail();
    }

  /* See if the intermediates are to be computed and set global variables
   * accordingly. */

  // this size estimate is only accurate if all centers are the same
  int size_inter_1 = cs1->nshell() * (sizeof(double*)+sizeof(int));
  if (storage - used_storage_ >= size_inter_1) {
      int_store1 = 1;
      used_storage_ += size_inter_1;
    }
  else {
    ExEnv::out0() << indent
         << "Int2eV3: not storing O(N) intemediates due to lack of memory"
         << endl;
    int_store1 = 0;
    }

  // this size estimate is only accurate if all centers are the same
  int size_inter_2 = cs1->nprimitive() * cs1->nprimitive() * (7*sizeof(double));
  if (storage - used_storage_ >= size_inter_2) {
      int_store2 = 1;
      used_storage_ += size_inter_2;
    }
  else {
    ExEnv::out0() << indent
         << "Int2eV3: not storing O(N^2) intermediates due to lack of memory"
         << endl;
    int_store2 = 0;
    }

  if (used_storage_ > storage || !int_store1 || !int_store2) {
    ExEnv::out0()
         << indent << "Int2eV3: wanted more storage than given" << endl
         << indent << "  given  storage = " << storage << endl
         << indent << "  build  storage = " << used_storage_build_ << endl
         << indent << "  shift  storage = " << used_storage_shift_ << endl
         << indent << "  used   storage = " << used_storage_ << endl
         << indent << "  O(N)   storage = " << size_inter_1
         <<           (int_store1?"":" (not used)") << endl
         << indent << "  O(N^2) storage = " << size_inter_2
         <<           (int_store2?"":" (not used)") << endl
         << endl;
    }

  int prim_inter_size = bs1_prim_offset_ + cs1->nprimitive();
  int shell_inter_size = bs1_shell_offset_ + cs1->nshell();
  if (bs2_prim_offset_ + (int_unit2?1:cs2->nprimitive()) > prim_inter_size) {
    prim_inter_size = bs2_prim_offset_ + (int_unit2?1:cs2->nprimitive());
    shell_inter_size = bs2_shell_offset_ + (int_unit2?1:cs2->nshell());
    }
  if (bs3_prim_offset_ + cs3->nprimitive() > prim_inter_size) {
    prim_inter_size = bs3_prim_offset_ + cs3->nprimitive();
    shell_inter_size = bs3_shell_offset_ + cs3->nshell();
    }
  if (bs4_prim_offset_ + (int_unit4?1:cs4->nprimitive()) > prim_inter_size) {
    prim_inter_size = bs4_prim_offset_ + (int_unit4?1:cs4->nprimitive());
    shell_inter_size = bs4_shell_offset_ + (int_unit4?1:cs4->nshell());
    }

  /* Allocate storage for the intermediates. */
  alloc_inter(prim_inter_size, shell_inter_size);

  /* Set up the one shell intermediates, block by block. */
  if (int_store1) {
    compute_shell_1(cs1, bs1_shell_offset_, bs1_prim_offset_);
    if (cs2.operator!=(cs1))
        compute_shell_1(cs2, bs2_shell_offset_, bs2_prim_offset_);
    if (cs3.operator!=(cs2) && cs3.operator!=(cs1))
        compute_shell_1(cs3, bs3_shell_offset_, bs3_prim_offset_);
    if (cs4.operator!=(cs3) && cs4.operator!=(cs2)&& cs4.operator!=(cs1))
        compute_shell_1(cs4, bs4_shell_offset_, bs4_prim_offset_);
    }

  /* Compute the two shell intermediates, block by block. */
  if (int_store2) {
    /* Compute the two primitive intermediates, block by block. */
    // Some unnecessary pairs of intermediates are avoided, but
    // some unnecessary pairs are still being computed.
    compute_prim_2(cs1,bs1_shell_offset_,bs1_prim_offset_,
                   cs1,bs1_shell_offset_,bs1_prim_offset_);
    if (cs2.operator!=(cs1)) {
      compute_prim_2(cs1,bs1_shell_offset_,bs1_prim_offset_,
                     cs2,bs2_shell_offset_,bs2_prim_offset_);
      compute_prim_2(cs2,bs2_shell_offset_,bs2_prim_offset_,
                     cs1,bs1_shell_offset_,bs1_prim_offset_);
      compute_prim_2(cs2,bs2_shell_offset_,bs2_prim_offset_,
                     cs2,bs2_shell_offset_,bs2_prim_offset_);
      }
    if (cs3.operator!=(cs2) && cs3.operator!=(cs1)) {
      compute_prim_2(cs1,bs1_shell_offset_,bs1_prim_offset_,
                     cs3,bs3_shell_offset_,bs3_prim_offset_);
      compute_prim_2(cs3,bs3_shell_offset_,bs3_prim_offset_,
                     cs1,bs1_shell_offset_,bs1_prim_offset_);
      compute_prim_2(cs2,bs2_shell_offset_,bs2_prim_offset_,
                     cs3,bs3_shell_offset_,bs3_prim_offset_);
      compute_prim_2(cs3,bs3_shell_offset_,bs3_prim_offset_,
                     cs2,bs2_shell_offset_,bs2_prim_offset_);
      compute_prim_2(cs3,bs3_shell_offset_,bs3_prim_offset_,
                     cs3,bs3_shell_offset_,bs3_prim_offset_);
      }
    if (cs4.operator!=(cs3) && cs4.operator!=(cs2) && cs4.operator!=(cs1)) {
      compute_prim_2(cs1,bs1_shell_offset_,bs1_prim_offset_,
                     cs4,bs4_shell_offset_,bs4_prim_offset_);
      compute_prim_2(cs4,bs4_shell_offset_,bs4_prim_offset_,
                     cs1,bs1_shell_offset_,bs1_prim_offset_);
      compute_prim_2(cs2,bs2_shell_offset_,bs2_prim_offset_,
                     cs4,bs4_shell_offset_,bs4_prim_offset_);
      compute_prim_2(cs4,bs4_shell_offset_,bs4_prim_offset_,
                     cs2,bs2_shell_offset_,bs2_prim_offset_);
      compute_prim_2(cs3,bs3_shell_offset_,bs3_prim_offset_,
                     cs4,bs4_shell_offset_,bs4_prim_offset_);
      compute_prim_2(cs4,bs4_shell_offset_,bs4_prim_offset_,
                     cs3,bs3_shell_offset_,bs3_prim_offset_);
      // cs4 cs4 terms are never needed since cs4 != cs3
      //compute_prim_2(cs4,bs4_shell_offset_,bs_prim_offset_,
      //               cs4,bs4_shell_offset_,bs_prim_offset_);
      }
    }

  return int_buffer;
  }

/* This is called when no more 2 electron integrals are needed.
 * It will free the intermediates. */
void
Int2eV3::int_done_erep()
{
  if (int_unit_shell) delete_int_unit_shell();
  if (int_derint_buffer) deallocate(int_derint_buffer);
  deallocate(int_buffer);
  if (int_store1) {
    delete[] int_shell_to_prim;
    }
  int_done_buildgc();
  int_done_shiftgc();
}

/* Allocates storage for the intermediates.  The arguments are the
 * total number of unique primitive and shells. */
void
Int2eV3::alloc_inter(int nprim,int nshell)
{
  if (int_store1) {
    int_shell_r.set_dim(nshell,3);
    int_shell_to_prim = new int[nshell];
    if (int_shell_to_prim == 0) {
      ExEnv::errn() << "problem allocating O(n) integral intermediates for";
      ExEnv::errn() << scprintf(" %d shells and %d primitives",nshell,nprim);
      ExEnv::errn() << endl;
      fail();
      }
    }
  if (int_store2) {
    int_prim_zeta.set_dim(nprim,nprim);
    int_prim_oo2zeta.set_dim(nprim,nprim);
    int_prim_k.set_dim(nprim,nprim);
    int_prim_p.set_dim(nprim,nprim,3);
    }
  }

void
Int2eV3::compute_shell_1(Ref<GaussianBasisSet> cs,
                         int shell_offset, int prim_offset)
{
  if (cs.null()) {
    for (int i=0; i<3; i++) {
      int_shell_r(shell_offset,i) = 0.0;
      }
    int_shell_to_prim[shell_offset] = prim_offset;
    return;
    }

  int i,j;
  int offset;
  int iprim;

  offset = shell_offset;
  iprim = prim_offset;
  for (i=0; i<cs->ncenter(); i++) {
    for (j=0; j<cs->nshell_on_center(i); j++) {

      /* The offset shell geometry vectors. */
      for (int xyz=0; xyz<3; xyz++) {
        int_shell_r(offset,xyz) = cs->molecule()->r(i,xyz);
        }

      /* The number of the first offset primitive in a offset shell. */
      int_shell_to_prim[offset] = iprim;

      offset++;
      iprim += cs->shell(i,j).nprimitive();
      }
    }
  }

/* The 2 primitive intermediates. */
void
Int2eV3::compute_prim_2(Ref<GaussianBasisSet> cs1,
                        int shell_offset1, int prim_offset1,
                        Ref<GaussianBasisSet> cs2,
                        int shell_offset2, int prim_offset2)
{
  int offset1, offset2;
  int i1,j1,k1,i2,j2,k2;
  GaussianShell *shell1,*shell2;
  int i;
  /* This is 2^(1/2) * pi^(5/4) */
  const double sqrt2pi54 = 5.9149671727956129;
  double AmB,AmB2;

  if (cs2.null() && !int_unit_shell) make_int_unit_shell();

  offset1 = prim_offset1;
  int cs1_ncenter = (cs1.null()?1:cs1->ncenter());
  for (i1=0; i1<cs1_ncenter; i1++) {
    int cs1_nshell_on_center = (cs1.null()?1:cs1->nshell_on_center(i1));
    for (j1=0; j1<cs1_nshell_on_center; j1++) {
      if (cs1.nonnull()) shell1 = &cs1->shell(i1,j1);
      else               shell1 = int_unit_shell;
      for (k1=0; k1<shell1->nprimitive(); k1++) {
        offset2 = prim_offset2;
        int cs2_ncenter = (cs2.null()?1:cs2->ncenter());
        for (i2=0; i2<cs2_ncenter; i2++) {
          int cs2_nshell_on_center = (cs2.null()?1:cs2->nshell_on_center(i2));
          for (j2=0; j2<cs2_nshell_on_center; j2++) {
            if (cs2.nonnull()) shell2 = &cs2->shell(i2,j2);
            else               shell2 = int_unit_shell;
            for (k2=0; k2<shell2->nprimitive(); k2++) {

              /* The zeta = alpha + beta intermediate. */
              int_prim_zeta(offset1,offset2) =
                shell1->exponent(k1) + shell2->exponent(k2);

              /* The 1/(2 zeta) intermediate times 2.0. */
              int_prim_oo2zeta(offset1,offset2) =
                1.0/int_prim_zeta(offset1,offset2);

              /* The p = (alpha A + beta B) / zeta */
              for (i=0; i<3; i++) {
                int_prim_p(offset1,offset2,i) =
                  (  shell1->exponent(k1) * (cs1.null()?0.0
                                             :cs1->molecule()->r(i1,i))
                   + shell2->exponent(k2) * (cs2.null()?0.0
                                             :cs2->molecule()->r(i2,i)))
                  *  int_prim_oo2zeta(offset1,offset2);
                }

              /* Compute AmB^2 */
              AmB2 = 0.0;
              for (i=0; i<3; i++) {
                AmB = (cs2.null()?0.0:cs2->molecule()->r(i2,i))
                    - (cs1.null()?0.0:cs1->molecule()->r(i1,i));
                AmB2 += AmB*AmB;
                }

              /* Compute the K intermediate. */
              int_prim_k(offset1,offset2) =
                   sqrt2pi54
                 * int_prim_oo2zeta(offset1,offset2)
                 * exp( -   shell1->exponent(k1) * shell2->exponent(k2)
                          * int_prim_oo2zeta(offset1,offset2)
                          * AmB2 );

              /* Finish the 1/(2 zeta) intermediate. */
              int_prim_oo2zeta(offset1,offset2) =
                0.5 * int_prim_oo2zeta(offset1,offset2);

              offset2++;
              }
            }
          }
        offset1++;
        }
      }
    }
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
