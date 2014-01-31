//
// comp1e.cc
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
#include <util/misc/math.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/tformv3.h>

using namespace std;
using namespace sc;

#define IN(i,j) ((i)==(j)?1:0)
#define SELECT(x1,x2,x3,s) (((s)==0)?x1:(((s)==1)?(x2):(x3)))

#define DEBUG_NUC_SHELL_DER 0
#define DEBUG_NUC_PRIM 0
#define DEBUG_NUC_SHELL 0

#define DEBUG_EFIELD_PRIM 0

/* ------------ Initialization of 1e routines. ------------------- */
/* This routine returns a buffer large enough to hold a shell doublet
 * of integrals (if order == 0) or derivative integrals (if order == 1).
 */
void
Int1eV3::int_initialize_1e(int flags, int order)
{
  int jmax1,jmax2,jmax;
  int scratchsize=0,nshell2;

  /* The efield routines look like derivatives.  The p_dot_nuclear_p
   * routines look like first derivatives in terms of the scratch buffer
   * size and second derivatives in terms of the primitive intermediates
   * needed.  Bump up order if it is zero or one to allow these integrals
   * to be computed.
   */
  int scratch_order = order;
  if (order <= 1) order = 2;
  if (scratch_order == 0) scratch_order = 1;

  jmax1 = bs1_->max_angular_momentum();
  jmax2 = bs2_->max_angular_momentum();
  jmax = jmax1 + jmax2;

  fjt_ = new FJT(jmax + 2*order);

  nshell2 = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell();

  if (scratch_order == 0) {
    init_order = 0;
    scratchsize = nshell2;
    }
  else if (scratch_order == 1) {
    init_order = 1;
    scratchsize = nshell2*3;
    }
  else {
    ExEnv::errn()
      << scprintf("int_initialize_1e: invalid scratch order: %d\n",
                  scratch_order);
    exit(1);
    }

  buff = new double[scratchsize];
  cartesianbuffer = new double[scratchsize];
  cartesianbuffer_scratch = new double[scratchsize];

  inter.set_dim(jmax1+order+1,jmax2+order+1,jmax+2*order+1);
  efield_inter.set_dim(jmax1+order+1,jmax2+order+1,jmax+2*order+1);
  int i,j,m;
  for (i=0; i<=jmax1+order; i++) {
    int sizei = INT_NCART_NN(i);
    for (j=0; j<=jmax2+order; j++) {
      int sizej = INT_NCART_NN(j);
      for (m=0; m<=jmax+2*order-i-j; m++) {
        inter(i,j,m) = new double[sizei*sizej];
        efield_inter(i,j,m) = new double[sizei*sizej*3];
        }
      for (; m<=jmax+2*order; m++) {
        inter(i,j,m) = 0;
        efield_inter(i,j,m) = 0;
        }
      }
    }

  }

void
Int1eV3::int_done_1e()
{
  init_order = -1;
  delete[] buff;
  delete[] cartesianbuffer;
  delete[] cartesianbuffer_scratch;
  buff = 0;
  cartesianbuffer = 0;
  inter.delete_data();
  efield_inter.delete_data();
}


/* --------------------------------------------------------------- */
/* ------------- Routines for the overlap matrix ----------------- */
/* --------------------------------------------------------------- */

/* This computes the overlap integrals between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::overlap(int ish, int jsh)
{
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int gc1,gc2;
  int index,index1,index2;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;
  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = comp_shell_overlap(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)

  transform_1e(integral_, cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the overlap ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::overlap_1der(int ish, int jsh,
                      int idercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_overlap: one electron routines are not init'ed\n");
    exit(1);
    }

  Ref<GaussianBasisSet> dercs;
  if (idercs == 0) dercs = bs1_;
  else dercs = bs2_;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

#if 0
  ExEnv::outn() << scprintf("zeroing %d*%d*3 elements of buff\n",ni,nj);
#endif
  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_overlap_1der(ish,jsh,dercs,centernum);
  }

/* This computes the overlap derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
void
Int1eV3::int_accum_shell_overlap_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs, int centernum)
{
  accum_shell_1der(buff,ish,jsh,dercs,centernum,&Int1eV3::comp_shell_overlap);
  }

/* Compute the overlap for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
double
Int1eV3::comp_shell_overlap(int gc1, int i1, int j1, int k1,
                            int gc2, int i2, int j2, int k2)
{
  double exp1,exp2;
  int i,j,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  if ((i1<0)||(j1<0)||(k1<0)||(i2<0)||(j2<0)||(k2<0)) return 0.0;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      exp1 = gshell1->exponent(i);
      exp2 = gshell2->exponent(j);
      oozeta = 1.0/(exp1 + exp2);
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(exp1 * A[xyz] + exp2 * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      ss =   pow(M_PI/(exp1+exp2),1.5)
           * exp(- oozeta * exp1 * exp2 * AmB2);
      tmp     =  gshell1->coefficient_unnorm(gc1,i)
               * gshell2->coefficient_unnorm(gc2,j)
               * comp_prim_overlap(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      result += tmp;
      }
    }

  return result;
  }

/* Compute the overlap between two primitive functions. */
#if 0
double
Int1eV3::int_prim_overlap(shell_t *pshell1, shell_t *pshell2,
                          double *pA, double *pB,
                          int prim1, int prim2,
                          int i1, int j1, int k1,
                          int i2, int j2, int k2)
{
  int xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;

  /* Compute the intermediates. */
  oozeta = 1.0/(gshell1->exponent(prim1) + gshell2->exponent(prim2));
  oo2zeta = 0.5*oozeta;
  AmB2 = 0.0;
  for (xyz=0; xyz<3; xyz++) {
    Pi = oozeta*(gshell1->exponent(prim1) * A[xyz]
                 + gshell2->exponent(prim2) * B[xyz]);
    PmA[xyz] = Pi - A[xyz];
    PmB[xyz] = Pi - B[xyz];
    AmB = A[xyz] - B[xyz];
    AmB2 += AmB*AmB;
    }
  ss =   pow(M_PI/(gshell1->exponent(prim1)
                                +gshell2->exponent(prim2)),1.5)
       * exp(- oozeta * gshell1->exponent(prim1)
             * gshell2->exponent(prim2) * AmB2);
  return comp_prim_overlap(i1,j1,k1,i2,j2,k2);
  }
#endif

double
Int1eV3::comp_prim_overlap(int i1, int j1, int k1,
                         int i2, int j2, int k2)
{
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_overlap(i1-1,j1,k1,i2,j2,k2);
    if (i1>1) result += oo2zeta*(i1-1) * comp_prim_overlap(i1-2,j1,k1,i2,j2,k2);
    if (i2>0) result += oo2zeta*i2 * comp_prim_overlap(i1-1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_overlap(i1,j1-1,k1,i2,j2,k2);
    if (j1>1) result += oo2zeta*(j1-1) * comp_prim_overlap(i1,j1-2,k1,i2,j2,k2);
    if (j2>0) result += oo2zeta*j2 * comp_prim_overlap(i1,j1-1,k1,i2,j2-1,k2);
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2);
    if (k1>1) result += oo2zeta*(k1-1) * comp_prim_overlap(i1,j1,k1-2,i2,j2,k2);
    if (k2>0) result += oo2zeta*k2 * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2-1);
    return result;
    }

  if (i2) {
    result = PmB[0] * comp_prim_overlap(i1,j1,k1,i2-1,j2,k2);
    if (i2>1) result += oo2zeta*(i2-1) * comp_prim_overlap(i1,j1,k1,i2-2,j2,k2);
    if (i1>0) result += oo2zeta*i1 * comp_prim_overlap(i1-1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_overlap(i1,j1,k1,i2,j2-1,k2);
    if (j2>1) result += oo2zeta*(j2-1) * comp_prim_overlap(i1,j1,k1,i2,j2-2,k2);
    if (j1>0) result += oo2zeta*j1 * comp_prim_overlap(i1,j1-1,k1,i2,j2-1,k2);
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_overlap(i1,j1,k1,i2,j2,k2-1);
    if (k2>1) result += oo2zeta*(k2-1) * comp_prim_overlap(i1,j1,k1,i2,j2,k2-2);
    if (k1>0) result += oo2zeta*k1 * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2-1);
    return result;
    }

  return ss;
  }

/* --------------------------------------------------------------- */
/* ------------- Routines for the kinetic energy ----------------- */
/* --------------------------------------------------------------- */

/* This computes the kinetic energy integrals between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::kinetic(int ish, int jsh)
{
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int cart1,cart2;
  int index;
  int gc1,gc2;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;
  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)

  transform_1e(integral_, cartesianbuffer, buff, gshell1, gshell2);
  }

void
Int1eV3::int_accum_shell_kinetic(int ish, int jsh)
{
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int cart1,cart2;
  int index;
  int gc1,gc2;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;

  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)
  accum_transform_1e(integral_, cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the kinetic energy derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
void
Int1eV3::int_accum_shell_kinetic_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs, int centernum)
{
  accum_shell_1der(buff,ish,jsh,dercs,centernum,&Int1eV3::comp_shell_kinetic);
  }

/* This computes the basis function part of 
 * generic derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 * The function used to compute the nonderivative integrals is shell_function.
 */
void
Int1eV3::accum_shell_1der(double *buff, int ish, int jsh,
                          Ref<GaussianBasisSet> dercs, int centernum,
                          double (Int1eV3::*shell_function)
                          (int,int,int,int,int,int,int,int))
{
  int i;
  int gc1,gc2;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index1,index2;
  double tmp[3];
  double *ctmp = cartesianbuffer;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      if ((bs1_==bs2_)&&(c1==c2)) {
        if (    three_center
             && !((bs1_==third_centers)&&(c1==third_centernum))
             && ((bs1_==dercs)&&(c1==centernum))) {
          for (i=0; i<3; i++) {
            /* Derivative wrt first shell. */
            exponent_weighted = 0;
            tmp[i] = 2.0 *
               (this->*shell_function)(gc1,i1+IN(i,0),j1+IN(i,1),k1+IN(i,2),gc2,i2,j2,k2);
            exponent_weighted = -1;
            if (SELECT(i1,j1,k1,i)) {
              tmp[i] -= SELECT(i1,j1,k1,i) *
                (this->*shell_function)(gc1,i1-IN(i,0),j1-IN(i,1),k1-IN(i,2),gc2,i2,j2,k2);
              }
            /* Derviative wrt second shell. */
            exponent_weighted = 1;
            tmp[i] += 2.0 *
               (this->*shell_function)(gc1,i1,j1,k1,gc2,i2+IN(i,0),j2+IN(i,1),k2+IN(i,2));
            exponent_weighted = -1;
            if (SELECT(i2,j2,k2,i)) {
              tmp[i] -= SELECT(i2,j2,k2,i) *
                (this->*shell_function)(gc1,i1,j1,k1,gc2,i2-IN(i,0),j2-IN(i,1),k2-IN(i,2));
              }
            }
	  }
        else {
          /* If there are two centers and they are the same, then we
           * use translational invariance to get a net contrib of 0.0 */
          for (i=0; i<3; i++) tmp[i] = 0.0;
          }
        }
      else if ((bs1_==dercs)&&(c1==centernum)) {
        for (i=0; i<3; i++) {
          exponent_weighted = 0;
          tmp[i] = 2.0 *
             (this->*shell_function)(gc1,i1+IN(i,0),j1+IN(i,1),k1+IN(i,2),gc2,i2,j2,k2);
          exponent_weighted = -1;
          if (SELECT(i1,j1,k1,i)) {
            tmp[i] -= SELECT(i1,j1,k1,i) *
              (this->*shell_function)(gc1,i1-IN(i,0),j1-IN(i,1),k1-IN(i,2),gc2,i2,j2,k2);
            }
          }
        }
      else if ((bs2_==dercs)&&(c2==centernum)) {
        for (i=0; i<3; i++) {
          exponent_weighted = 1;
          tmp[i] = 2.0 *
             (this->*shell_function)(gc1,i1,j1,k1,gc2,i2+IN(i,0),j2+IN(i,1),k2+IN(i,2));
          exponent_weighted = -1;
          if (SELECT(i2,j2,k2,i)) {
            tmp[i] -= SELECT(i2,j2,k2,i) *
              (this->*shell_function)(gc1,i1,j1,k1,gc2,i2-IN(i,0),j2-IN(i,1),k2-IN(i,2));
            }
          }
        }
      else {
        for (i=0; i<3; i++) tmp[i] = 0.0;
        }

      if (scale_shell_result) {
        for (i=0; i<3; i++) tmp[i] *= result_scale_factor;
        }

      for (i=0; i<3; i++) ctmp[i] = tmp[i];

      /* Increment the pointer to the xyz for the next atom. */
      ctmp += 3;
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)

  accum_transform_1e_xyz(integral_,
                         cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the basis function part of 
 * generic derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 * The function used to compute the nonderivative integrals is shell_function.
 */
void
Int1eV3::accum_shell_block_1der(double *buff, int ish, int jsh,
                                Ref<GaussianBasisSet> dercs, int centernum,
                                void (Int1eV3::*shell_block_function)
                                  (int gc1, int a, int gc2, int b,
                                   int gcsize2, int gcoff1, int gcoff2,
                                   double coef, double *buffer))
{
  int i;
  int gc1,gc2;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index1,index2;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);

  int docenter1=0, docenter2=0;
  int equiv12 = (bs1_==bs2_)&&(c1==c2);
  int der1 = (bs1_==dercs)&&(c1==centernum);
  int der2 = (bs2_==dercs)&&(c2==centernum);
  if (!equiv12) {
    docenter1 = der1;
    docenter2 = der2;
    }
  else if (three_center) {
    int equiv123 = (bs1_==third_centers)&&(c1==third_centernum);
    if (!equiv123) {
      docenter1 = der1;
      docenter2 = der2;
      }
    }

  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  int gcsize1 = gshell1->ncartesian();
  int gcsize2 = gshell2->ncartesian();
  memset(cartesianbuffer,0,sizeof(double)*gcsize1*gcsize2*3);

  if (!docenter1 && !docenter2) return;

  double coef;
  if (scale_shell_result) {
    coef = result_scale_factor;
    }
  else coef = 1.0;

  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  int gcoff1 = 0;
  for (gc1=0; gc1<gshell1->ncontraction(); gc1++) {
    int a = gshell1->am(gc1);
    int sizea = INT_NCART_NN(a);
    int sizeap1 = INT_NCART_NN(a+1);
    int sizeam1 = INT_NCART(a-1);
    int gcoff2 = 0;
    for (gc2=0; gc2<gshell2->ncontraction(); gc2++) {
      int b = gshell2->am(gc2);
      int sizeb = INT_NCART_NN(b);
      int sizebp1 = INT_NCART_NN(b+1);
      int sizebm1 = INT_NCART(b-1);
      /* Derivative wrt first shell. */
      if (docenter1) {
        exponent_weighted = 0;
        memset(cartesianbuffer_scratch,0,sizeof(double)*sizeap1*sizeb);
        (this->*shell_block_function)(gc1, a+1, gc2, b,
                                      sizeb, 0, 0,
                                      coef, cartesianbuffer_scratch);
        index1=0;
        FOR_CART(i1,j1,k1,a) {
          index2=0;
          FOR_CART(i2,j2,k2,b) {
            double *ctmp = &cartesianbuffer[((index1+gcoff1)*gcsize2
                                             +index2+gcoff2)*3];
            for (i=0; i<3; i++) {
              int ind = INT_CARTINDEX(a+1,i1+IN(i,0),j1+IN(i,1));
              *ctmp++ += 2.0*cartesianbuffer_scratch[ind*sizeb+index2];
              }
            index2++;
            } END_FOR_CART;
          index1++;
          } END_FOR_CART;

        if (a) {
          exponent_weighted = -1;
          memset(cartesianbuffer_scratch,0,sizeof(double)*sizeam1*sizeb);
          (this->*shell_block_function)(gc1, a-1, gc2, b,
                                        sizeb, 0, 0,
                                        coef, cartesianbuffer_scratch);
          index1=0;
          FOR_CART(i1,j1,k1,a) {
            index2=0;
            FOR_CART(i2,j2,k2,b) {
              double *ctmp = &cartesianbuffer[((index1+gcoff1)*gcsize2
                                               +index2+gcoff2)*3];
              for (i=0; i<3; i++) {
                int sel = SELECT(i1,j1,k1,i);
                if (sel) {
                  int ind = INT_CARTINDEX(a-1,i1-IN(i,0),j1-IN(i,1));
                  ctmp[i] -= sel * cartesianbuffer_scratch[ind*sizeb+index2];
                  }
                }
              index2++;
              } END_FOR_CART;
            index1++;
            } END_FOR_CART;
          }
        }
      if (docenter2) {
        /* Derviative wrt second shell. */
        exponent_weighted = 1;
        memset(cartesianbuffer_scratch,0,sizeof(double)*sizea*sizebp1);
        (this->*shell_block_function)(gc1, a, gc2, b+1,
                                      sizebp1, 0, 0,
                                      coef, cartesianbuffer_scratch);
        index1=0;
        FOR_CART(i1,j1,k1,a) {
          index2=0;
          FOR_CART(i2,j2,k2,b) {
            double *ctmp = &cartesianbuffer[((index1+gcoff1)*gcsize2
                                             +index2+gcoff2)*3];
            for (i=0; i<3; i++) {
              int ind = INT_CARTINDEX(b+1,i2+IN(i,0),j2+IN(i,1));
              *ctmp++ += 2.0*cartesianbuffer_scratch[index1*sizebp1+ind];
              }
            index2++;
            } END_FOR_CART;
          index1++;
          } END_FOR_CART;

        if (b) {
          exponent_weighted = -1;
          memset(cartesianbuffer_scratch,0,sizeof(double)*sizea*sizebm1);
          (this->*shell_block_function)(gc1, a, gc2, b-1,
                                        sizebm1, 0, 0,
                                        coef, cartesianbuffer_scratch);
          index1=0;
          FOR_CART(i1,j1,k1,a) {
            index2=0;
            FOR_CART(i2,j2,k2,b) {
              double *ctmp = &cartesianbuffer[((index1+gcoff1)*gcsize2
                                               +index2+gcoff2)*3];
              for (i=0; i<3; i++) {
                int sel = SELECT(i2,j2,k2,i);
                if (sel) {
                  int ind = INT_CARTINDEX(b-1,i2-IN(i,0),j2-IN(i,1));
                  ctmp[i] -= sel * cartesianbuffer_scratch[index1*sizebm1+ind];
                  }
                }
              index2++;
              } END_FOR_CART;
            index1++;
            } END_FOR_CART;
          }
        }
      gcoff2 += INT_NCART_NN(b);
      }
    gcoff1 += INT_NCART_NN(a);
    }

  accum_transform_1e_xyz(integral_,
                         cartesianbuffer, buff, gshell1, gshell2);

#if DEBUG_NUC_SHELL_DER
  double *fastbuff = cartesianbuffer;
  cartesianbuffer = new double[gcsize1*gcsize2*3];
  double *junkbuff = new double[gcsize1*gcsize2*3];
  memset(junkbuff,0,sizeof(double)*gcsize1*gcsize2*3);
  accum_shell_1der(junkbuff,ish,jsh,dercs,centernum,
                   &Int1eV3::comp_shell_nuclear);
  delete[] junkbuff;
  int index = 0;
  for (i=0; i<gcsize1; i++) {
      for (int j=0; j<gcsize2; j++, index++) {
          ExEnv::outn() << scprintf(" (%d %d): % 18.15f % 18.15f",
                           i,j,cartesianbuffer[index],fastbuff[index]);
          if (fabs(cartesianbuffer[index]-fastbuff[index])>1.0e-12) {
              ExEnv::outn() << " **";
            }
          ExEnv::outn() << endl;
        }
    }
  delete[] cartesianbuffer;
  cartesianbuffer = fastbuff;
#endif
  }

/* Compute the kinetic energy for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
double
Int1eV3::comp_shell_kinetic(int gc1, int i1, int j1, int k1,
                          int gc2, int i2, int j2, int k2)
{
  int i,j,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      oo2zeta_a = 0.5/gshell1->exponent(i);
      oo2zeta_b = 0.5/gshell2->exponent(j);
      oozeta = 1.0/(gshell1->exponent(i) + gshell2->exponent(j));
      oo2zeta = 0.5*oozeta;
      xi = oozeta * gshell1->exponent(i) * gshell2->exponent(j);
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz]
                     + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      /* The s integral kinetic energy. */
      ss =   pow(M_PI/(gshell1->exponent(i)
                                    +gshell2->exponent(j)),1.5)
           * exp(- xi * AmB2);
      sTs =  ss
           * xi
           * (3.0 - 2.0 * xi * AmB2);
      tmp     =  gshell1->coefficient_unnorm(gc1,i)
               * gshell2->coefficient_unnorm(gc2,j)
               * comp_prim_kinetic(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= gshell1->exponent(i);
      else if (exponent_weighted == 1) tmp *= gshell2->exponent(j);
      result += tmp;
      }
    }

  return result;
  }

double
Int1eV3::comp_prim_kinetic(int i1, int j1, int k1,
                         int i2, int j2, int k2)
{
  double tmp;
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_kinetic(i1-1,j1,k1,i2,j2,k2);
    if (i1>1) result += oo2zeta*(i1-1)*comp_prim_kinetic(i1-2,j1,k1,i2,j2,k2);
    if (i2) result += oo2zeta*i2*comp_prim_kinetic(i1-1,j1,k1,i2-1,j2,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (i1>1) tmp -= oo2zeta_a*(i1-1)*comp_prim_overlap(i1-2,j1,k1,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_kinetic(i1,j1-1,k1,i2,j2,k2);
    if (j1>1) result += oo2zeta*(j1-1)*comp_prim_kinetic(i1,j1-2,k1,i2,j2,k2);
    if (j2) result += oo2zeta*j2*comp_prim_kinetic(i1,j1-1,k1,i2,j2-1,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (j1>1) tmp -= oo2zeta_a*(j1-1)*comp_prim_overlap(i1,j1-2,k1,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2);
    if (k1>1) result += oo2zeta*(k1-1)*comp_prim_kinetic(i1,j1,k1-2,i2,j2,k2);
    if (k2) result += oo2zeta*k2*comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2-1);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (k1>1) tmp -= oo2zeta_a*(k1-1)*comp_prim_overlap(i1,j1,k1-2,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (i2) {
    result = PmB[0] * comp_prim_kinetic(i1,j1,k1,i2-1,j2,k2);
    if (i2>1) result += oo2zeta*(i2-1)*comp_prim_kinetic(i1,j1,k1,i2-2,j2,k2);
    if (i1) result += oo2zeta*i1*comp_prim_kinetic(i1-1,j1,k1,i2-1,j2,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (i2>1) tmp -= oo2zeta_b*(i2-1)*comp_prim_overlap(i1,j1,k1,i2-2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_kinetic(i1,j1,k1,i2,j2-1,k2);
    if (j2>1) result += oo2zeta*(j2-1)*comp_prim_kinetic(i1,j1,k1,i2,j2-2,k2);
    if (j1) result += oo2zeta*i1*comp_prim_kinetic(i1,j1-1,k1,i2,j2-1,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (j2>1) tmp -= oo2zeta_b*(j2-1)*comp_prim_overlap(i1,j1,k1,i2,j2-2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_kinetic(i1,j1,k1,i2,j2,k2-1);
    if (k2>1) result += oo2zeta*(k2-1)*comp_prim_kinetic(i1,j1,k1,i2,j2,k2-2);
    if (k1) result += oo2zeta*i1*comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2-1);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (k2>1) tmp -= oo2zeta_b*(k2-1)*comp_prim_overlap(i1,j1,k1,i2,j2,k2-2);
    result += 2.0 * xi * tmp;
    return result;
    }

  return sTs;
  }

/* --------------------------------------------------------------- */
/* ------------- Routines for the nuclear attraction ------------- */
/* --------------------------------------------------------------- */

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
void
Int1eV3::int_accum_shell_nuclear_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs, int centernum)
{
  int_accum_shell_nuclear_hf_1der(ish,jsh,dercs,centernum);
  int_accum_shell_nuclear_nonhf_1der(ish,jsh,dercs,centernum);
  }

/* A correction to the Hellman-Feynman part is computed which
 * is not included in the original HF routine.  This is only needed
 * if the real Hellman-Feynman forces are desired, because the sum
 * of the hf_1der and nonhf_1der forces are still correct.
 */
void
Int1eV3::int_accum_shell_nuclear_hfc_1der(int ish, int jsh,
                                          Ref<GaussianBasisSet> dercs,
                                          int centernum)
{
  /* If both ish and jsh are not on the der center,
   * then there's no correction. */
  if (!(  (bs1_==dercs)
        &&(bs2_==dercs)
        &&(bs1_->shell_to_center(ish)==centernum)
        &&(bs2_->shell_to_center(jsh)==centernum))) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative for three equal
   * centers. */
  scale_shell_result = 1;
  result_scale_factor = -bs1_->molecule()->charge(centernum);
  for (int xyz=0; xyz<3; xyz++) {
      C[xyz] = bs1_->r(centernum,xyz);
    }
  accum_shell_efield(buff,ish,jsh);
  scale_shell_result = 0;

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the Hellman-Feynman part is computed.
 */
void
Int1eV3::int_accum_shell_nuclear_hf_1der(int ish, int jsh,
                                         Ref<GaussianBasisSet> dercs,
                                         int centernum)
{

  /* If both ish and jsh are on the der center, then the contrib is zero. */
  if (  (bs1_==dercs)
      &&(bs2_==dercs)
      &&(bs1_->shell_to_center(ish)==centernum)
      &&(bs2_->shell_to_center(jsh)==centernum)) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative. */
  if (bs1_ == dercs) {
    scale_shell_result = 1;
    result_scale_factor= -bs1_->molecule()->charge(centernum);
    for (int xyz=0; xyz<3; xyz++) {
        C[xyz] = bs1_->r(centernum,xyz);
      }
    //accum_shell_efield(buff,ish,jsh);
    accum_shell_block_efield(buff,ish,jsh);
    scale_shell_result = 0;
    }
  else if (bs2_ == dercs) {
    scale_shell_result = 1;
    result_scale_factor= -bs2_->molecule()->charge(centernum);
    for (int xyz=0; xyz<3; xyz++) {
        C[xyz] = bs2_->r(centernum,xyz);
      }
    //accum_shell_efield(buff,ish,jsh);
    accum_shell_block_efield(buff,ish,jsh);
    scale_shell_result = 0;
    }

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the non Hellman-Feynman part is computed.
 */
void
Int1eV3::int_accum_shell_nuclear_nonhf_1der(int ish, int jsh,
                                            Ref<GaussianBasisSet> dercs,
                                            int centernum)
{
  int i;

  /* Get the basis function part of the nuclear derivative. */
  three_center = 1;
  third_centers = bs1_;
  for (i=0; i<bs1_->ncenter(); i++) {
    third_centernum = i;
    for (int xyz=0; xyz<3; xyz++) {
        C[xyz] = bs1_->r(i,xyz);
      }
    scale_shell_result = 1;
    result_scale_factor = -bs1_->molecule()->charge(i);
    //accum_shell_1der(buff,ish,jsh,dercs,centernum,
    //                 &Int1eV3::comp_shell_nuclear);
    accum_shell_block_1der(buff,ish,jsh,dercs,centernum,
                           &Int1eV3::comp_shell_block_nuclear);
    scale_shell_result = 0;
    }
  if (bs2_!=bs1_) {
    third_centers = bs2_;
    for (i=0; i<bs2_->ncenter(); i++) {
      third_centernum = i;
      for (int xyz=0; xyz<3; xyz++) {
          C[xyz] = bs2_->r(i,xyz);
        }
      scale_shell_result = 1;
      result_scale_factor = -bs2_->molecule()->charge(i);
      //accum_shell_1der(buff,ish,jsh,dercs,centernum,
      //                 &Int1eV3::comp_shell_nuclear);
      accum_shell_block_1der(buff,ish,jsh,dercs,centernum,
                             &Int1eV3::comp_shell_block_nuclear);
      scale_shell_result = 0;
      }
    }
  three_center = 0;

  }

/* This computes the efield integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.
 */
void
Int1eV3::int_accum_shell_efield(int ish, int jsh,
                                double *position)
{
  scale_shell_result = 0;
  for (int xyz=0; xyz<3; xyz++) {
      C[xyz] = position[xyz];
    }
  accum_shell_efield(buff,ish,jsh);
}

/* This computes the efield integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The globals scale_shell_result, result_scale_factor,
 * and C must be set before this is called.
 */
void
Int1eV3::accum_shell_efield(double *buff, int ish, int jsh)
{
  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    ExEnv::errn() << scprintf("accum_shell_efield: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2);
      if (scale_shell_result) {
        for (i=0; i<3; i++) efield[i] *= result_scale_factor;
        }
      for (i=0; i<3; i++) tmp[i] = efield[i];
      tmp += 3;
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)

  accum_transform_1e_xyz(integral_,
                         cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the efield integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The globals scale_shell_result, result_scale_factor,
 * and C must be set before this is called.
 */
void
Int1eV3::accum_shell_block_efield(double *buff, int ish, int jsh)
{
  int c1,c2;
  int gc1,gc2;

  if (!(init_order >= 1)) {
    ExEnv::errn() << scprintf("accum_shell_efield: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  double coef;
  if (scale_shell_result) coef = result_scale_factor;
  else coef = 1.0;

  int gcsize1 = gshell1->ncartesian();
  int gcsize2 = gshell2->ncartesian();
  memset(cartesianbuffer,0,sizeof(double)*gcsize1*gcsize2*3);
  int gcoff1 = 0;
  for (gc1=0; gc1<gshell1->ncontraction(); gc1++) {
    int a = gshell1->am(gc1);
    int sizea = INT_NCART_NN(a);
    int gcoff2 = 0;
    for (gc2=0; gc2<gshell2->ncontraction(); gc2++) {
      int b = gshell2->am(gc2);
      int sizeb = INT_NCART_NN(b);
      comp_shell_block_efield(gc1,a,gc2,b,
                              gcsize2, gcoff1, gcoff2,
                              coef, cartesianbuffer);
      gcoff2 += sizeb;
      }
    gcoff1 += sizea;
    }

  accum_transform_1e_xyz(integral_,
                         cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the efield integrals between functions in two shells.
 * The result is placed in the buffer in the form bf1 x y z, bf2
 * x y z, etc.
 */
void
Int1eV3::efield(int ish, int jsh, const double *position)
{
  scale_shell_result = 0;
  int xyz;

  for (xyz=0; xyz<3; xyz++) {
      C[xyz] = position[xyz];
    }

  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    ExEnv::errn() << scprintf("Int1eV3::efield one electron routines are not ready\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2);
      if (scale_shell_result) {
        for (i=0; i<3; i++) efield[i] *= result_scale_factor;
        }
      for (i=0; i<3; i++) tmp[i] = efield[i];
      tmp += 3;
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)

  transform_1e_xyz(integral_,
                   cartesianbuffer, buff, gshell1, gshell2);
}

/* This slowly computes the nuc rep energy integrals between functions in
 * two shells.  The result is placed in the buffer.  */
void
Int1eV3::nuclear_slow(int ish, int jsh)
{
  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;

  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = 0.0;
      /* Loop thru the centers on bs1_. */
      for (i=0; i<bs1_->ncenter(); i++) {
        for (int xyz=0; xyz<3; xyz++) {
          C[xyz] = bs1_->r(i,xyz);
          }
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * bs1_->molecule()->charge(i);
        }
      /* Loop thru the centers on bs2_ if necessary. */
      if (bs2_ != bs1_) {
        for (i=0; i<bs2_->ncenter(); i++) {
          for (int xyz=0; xyz<3; xyz++) {
            C[xyz] = bs2_->r(i,xyz);
            }
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                         * bs2_->molecule()->charge(i);
          }
        }
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)

  transform_1e(integral_,
               cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the nuc rep energy integrals between functions in two
 * shells.  The result is placed in the buffer.  */
void
Int1eV3::nuclear(int ish, int jsh)
{
  int i;
  int c1,c2;
  int gc1,gc2;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  int ni = gshell1->ncartesian();
  int nj = gshell2->ncartesian();
  memset(cartesianbuffer,0,sizeof(double)*ni*nj);

  int offi = 0;
  for (gc1=0; gc1<gshell1->ncontraction(); gc1++) {
    int a = gshell1->am(gc1);
    int offj = 0;
    for (gc2=0; gc2<gshell2->ncontraction(); gc2++) {
      int b = gshell2->am(gc2);
      /* Loop thru the centers on bs1_. */
      for (i=0; i<bs1_->ncenter(); i++) {
        double charge = bs1_->molecule()->charge(i);
        for (int xyz=0; xyz<3; xyz++) C[xyz] = bs1_->r(i,xyz);
        comp_shell_block_nuclear(gc1, a, gc2, b,
                                 nj, offi, offj,
                                 -charge, cartesianbuffer);
        }
      /* Loop thru the centers on bs2_ if necessary. */
      if (bs2_ != bs1_) {
        for (i=0; i<bs2_->ncenter(); i++) {
          double charge = bs2_->molecule()->charge(i);
          for (int xyz=0; xyz<3; xyz++) C[xyz] = bs2_->r(i,xyz);
          comp_shell_block_nuclear(gc1, a, gc2, b,
                                   nj, offi, offj,
                                   -charge, cartesianbuffer);
          }
        }
      offj += INT_NCART_NN(b);
      }
    offi += INT_NCART_NN(a);
    }

#if DEBUG_NUC_SHELL
  double *fastbuf = new double[ni*nj];
  memcpy(fastbuf,cartesianbuffer,sizeof(double)*ni*nj);
  nuclear_slow(ish,jsh);

  index = 0;
  int cart1, cart2;
  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1) {
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2) {
      double fast = fastbuf[index];
      double slow = cartesianbuffer[index];
      if (fabs(fast-slow)>1.0e-12) {
        ExEnv::outn() << scprintf("NUC SHELL FINAL: %d (%d %d %d) %d (%d %d %d): ",
                         gc1, i1,j1,k1, gc2, i2,j2,k2)
             << scprintf(" % 20.15f % 20.15f",
                         fast, slow)
             << endl;
        }
      index++;
      } END_FOR_GCCART_GS(cart2);
    } END_FOR_GCCART_GS(cart1);

  memcpy(cartesianbuffer,fastbuf,sizeof(double)*ni*nj);
  delete[] fastbuf;
#endif

  transform_1e(integral_,
               cartesianbuffer, buff, gshell1, gshell2);
  }

void
Int1eV3::p_dot_nuclear_p(int ish, int jsh)
{
  int i;
  int c1,c2;
  int gc1,gc2;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("Int1eV3::p_dot_nuclear_p: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  int ni = gshell1->ncartesian();
  int nj = gshell2->ncartesian();
  memset(cartesianbuffer,0,sizeof(double)*ni*nj);

  int offi = 0;
  for (gc1=0; gc1<gshell1->ncontraction(); gc1++) {
    int a = gshell1->am(gc1);
    int offj = 0;
    for (gc2=0; gc2<gshell2->ncontraction(); gc2++) {
      int b = gshell2->am(gc2);
      /* Loop thru the centers on bs1_. */
      for (i=0; i<bs1_->ncenter(); i++) {
        double charge = bs1_->molecule()->charge(i);
        for (int xyz=0; xyz<3; xyz++) C[xyz] = bs1_->r(i,xyz);
        comp_shell_block_p_dot_nuclear_p(gc1, a, gc2, b,
                                         nj, offi, offj,
                                         -charge, cartesianbuffer);
        }
      /* Loop thru the centers on bs2_ if necessary. */
      if (bs2_ != bs1_) {
        for (i=0; i<bs2_->ncenter(); i++) {
          double charge = bs2_->molecule()->charge(i);
          for (int xyz=0; xyz<3; xyz++) C[xyz] = bs2_->r(i,xyz);
          comp_shell_block_p_dot_nuclear_p(gc1, a, gc2, b,
                                           nj, offi, offj,
                                           -charge, cartesianbuffer);
          }
        }
      offj += INT_NCART_NN(b);
      }
    offi += INT_NCART_NN(a);
    }

  transform_1e(integral_,
               cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
void
Int1eV3::int_accum_shell_point_charge(int ish, int jsh,
                                      int ncharge, const double* charge,
                                      const double*const* position)
{
  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double tmp;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_pointcharge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;

  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      /* Loop thru the point charges. */
      tmp = 0.0;
      for (i=0; i<ncharge; i++) {
          for (int xyz=0; xyz<3; xyz++) {
              C[xyz] = position[i][xyz];
            }
        tmp -=  comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * charge[i];
        }
      cartesianbuffer[index] = tmp;
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)

  accum_transform_1e(integral_,
                     cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
void
Int1eV3::point_charge(int ish, int jsh,
                      int ncharge,
                      const double* charge, const double*const* position)
{
  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("Int1eV3::point_charge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;

  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = 0.0;
      /* Loop thru the point charges. */
      for (i=0; i<ncharge; i++) {
        for (int xyz=0; xyz<3; xyz++) {
          C[xyz] = position[i][xyz];
          }
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * charge[i];
        }
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)

  transform_1e(integral_,
               cartesianbuffer, buff, gshell1, gshell2);
  }


/* This computes the 1e Hamiltonian integrals between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::hcore(int ish, int jsh)
{
  int i;
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int index;
  int cart1,cart2;
  int gc1,gc2;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  index = 0;
  FOR_GCCART_GS(gc1,cart1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,cart2,i2,j2,k2,gshell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      /* Loop thru the centers on bs1_. */
      for (i=0; i<bs1_->ncenter(); i++) {
        for (int xyz=0; xyz<3; xyz++) {
          C[xyz] = bs1_->r(i,xyz);
          }
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * bs1_->molecule()->charge(i);
        }
      /* Loop thru the centers on bs2_ if necessary. */
      if (bs2_ != bs1_) {
        for (i=0; i<bs2_->ncenter(); i++) {
          for (int xyz=0; xyz<3; xyz++) {
            C[xyz] = bs2_->r(i,xyz);
            }
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * bs2_->molecule()->charge(i);
          }
        }
      index++;
      END_FOR_GCCART_GS(cart2)
    END_FOR_GCCART_GS(cart1)

  transform_1e(integral_,
               cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the 1e Hamiltonian deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::hcore_1der(int ish, int jsh,
                    int idercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  Ref<GaussianBasisSet> dercs;
  if (idercs == 0) dercs = bs1_;
  else dercs = bs2_;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_1der(ish,jsh,dercs,centernum);
  int_accum_shell_kinetic_1der(ish,jsh,dercs,centernum);
  }

/* This computes the kinetic deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::kinetic_1der(int ish, int jsh,
                      int idercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_kinetic: one electron routines are not init'ed\n");
    exit(1);
    }

  Ref<GaussianBasisSet> dercs;
  if (idercs == 0) dercs = bs1_;
  else dercs = bs2_;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_kinetic_1der(ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::nuclear_1der(int ish, int jsh, int idercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  Ref<GaussianBasisSet> dercs;
  if (idercs == 0) dercs = bs1_;
  else dercs = bs2_;

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_1der(ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * Only the Hellman-Feynman part is computed.
 * The result is placed in the buffer.
 */
void
Int1eV3::int_shell_nuclear_hf_1der(int ish, int jsh,
                                   Ref<GaussianBasisSet> dercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_nuclear_hf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_hf_1der(ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * Only the non Hellman-Feynman part is computed.
 * The result is placed in the buffer.
 */
void
Int1eV3::int_shell_nuclear_nonhf_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs, int centernum)
{
  int i;
  int c1,c2;
  int ni,nj;

  if (!(init_order >= 0)) {
    ExEnv::errn() << scprintf("int_shell_nuclear_nonhf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (int xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);

  ni = gshell1->nfunction();
  nj = gshell2->nfunction();

#if 0
  ExEnv::outn() << scprintf("int_shell_nuclear_nonhf_1der: zeroing %d doubles in buff\n",ni*nj*3);
#endif
  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_nonhf_1der(ish,jsh,dercs,centernum);
  }

/* Compute the nuclear attraction for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
double
Int1eV3::comp_shell_nuclear(int gc1, int i1, int j1, int k1,
                          int gc2, int i2, int j2, int k2)
{
  int i,j,k,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  int am;
  double tmp;

  am = i1+j1+k1+i2+j2+k2;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      zeta = gshell1->exponent(i) + gshell2->exponent(j);
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz]
                     + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * M_PI/(gshell1->exponent(i)
                                           +gshell2->exponent(j))
           * exp(- oozeta * gshell1->exponent(i)
                 * gshell2->exponent(j) * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am; k++) {
        fjttable_[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      tmp     =  gshell1->coefficient_unnorm(gc1,i)
               * gshell2->coefficient_unnorm(gc2,j)
               * comp_prim_nuclear(i1,j1,k1,i2,j2,k2,0);

      if (exponent_weighted == 0) tmp *= gshell1->exponent(i);
      else if (exponent_weighted == 1) tmp *= gshell2->exponent(j);

      result += tmp;
      }
    }

  /* printf("comp_shell_nuclear(%d,%d,%d,%d,%d,%d): result = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,result);
   */
  return result;
  }

/* Compute the nuclear attraction for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
void
Int1eV3::comp_shell_block_nuclear(int gc1, int a, int gc2, int b,
                                  int gcsize2, int gcoff1, int gcoff2,
                                  double coef, double *buffer)
{
  int i,j,k,xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  double tmp;
  int am = a + b;
  int size1 = INT_NCART_NN(a);
  int size2 = INT_NCART_NN(b);

#if DEBUG_NUC_SHELL
  double *shellints = new double[size1*size2];
  memset(shellints,0,sizeof(double)*size1*size2);
#endif

  /* Loop over the primitives in the shells. */
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      zeta = gshell1->exponent(i) + gshell2->exponent(j);
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz]
                     + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * M_PI/(gshell1->exponent(i)
                                           +gshell2->exponent(j))
           * exp(- oozeta * gshell1->exponent(i)
                 * gshell2->exponent(j) * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am; k++) {
        fjttable_[k] *= auxcoef;
        }

      /* Compute the primitive nuclear attraction integral. */
      comp_prim_block_nuclear(a,b);

      tmp =  gshell1->coefficient_unnorm(gc1,i)
             * gshell2->coefficient_unnorm(gc2,j)
             * coef;

      if (exponent_weighted == 0) tmp *= gshell1->exponent(i);
      else if (exponent_weighted == 1) tmp *= gshell2->exponent(j);

#if DEBUG_NUC_SHELL
      double *tprimbuffer = inter(a,b,0);
      double *tmpshellints = shellints;
      for (int ip=0; ip<size1; ip++) {
        for (int jp=0; jp<size2; jp++) {
          *tmpshellints++ += tmp * *tprimbuffer++ / coef;
          }
        }
#endif
      double *primbuffer = inter(a,b,0);
      for (int ip=0; ip<size1; ip++) {
        for (int jp=0; jp<size2; jp++) {
          //ExEnv::outn() << scprintf("buffer[%d] += %18.15f",
          //                 (ip+gcoff1)*gcsize2+jp+gcoff2,
          //                 tmp * *primbuffer)
          //     << endl;
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2] += tmp * *primbuffer++;
          }
        }
      }
    }

#if DEBUG_NUC_SHELL
#  if DEBUG_NUC_SHELL > 1
  ExEnv::outn() << scprintf("GC = (%d %d), A = %d, B = %d", gc1, gc2, a, b)
       << endl;
#  endif
  int i1,j1,k1;
  int i2,j2,k2;
  int ip = 0;
  double *tmpshellints = shellints;
  FOR_CART(i1,j1,k1,a) {
    int jp = 0;
    FOR_CART(i2,j2,k2,b) {
      double fast = *tmpshellints++;
      double slow = comp_shell_nuclear(gc1, i1, j1, k1,
                                       gc2, i2, j2, k2);
      int bad = fabs(fast-slow)>1.0e-12;
      if (DEBUG_NUC_SHELL > 1 || bad) {
        ExEnv::outn() << scprintf("NUC SHELL: (%d %d %d) (%d %d %d): ",
                         i1,j1,k1, i2,j2,k2)
             << scprintf(" % 20.15f % 20.15f",
                         fast, slow);
        }
      if (bad) {
        ExEnv::outn() << " ****" << endl;
        }
      else if (DEBUG_NUC_SHELL > 1) {
        ExEnv::outn() << endl;
        }
      jp++;
      } END_FOR_CART;
    ip++;
    } END_FOR_CART;
  delete[] shellints;
#endif
  }

void
Int1eV3::comp_shell_block_p_dot_nuclear_p(int gc1, int a, int gc2, int b,
                                          int gcsize2, int gcoff1, int gcoff2,
                                          double coef, double *buffer)
{
  int i,j,k,xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  double tmp;
  int am = a + b;
  int size1 = INT_NCART_NN(a);
  int size2 = INT_NCART_NN(b);

  /* Loop over the primitives in the shells. */
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      double alpha2 = gshell2->exponent(j);

      /* Compute the intermediates. */
      zeta = gshell1->exponent(i) + gshell2->exponent(j);
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz]
                     + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

// define this to compute only kinetic energy like integrals
// the computed integrals will be - 2 * T
#define PVP_KINETIC_ONLY 0
#define PVP_DEBUG 0

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * M_PI/(gshell1->exponent(i)
                                           +gshell2->exponent(j))
           * exp(- oozeta * gshell1->exponent(i)
                 * gshell2->exponent(j) * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am+2,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am+2; k++) {
        fjttable_[k] *= auxcoef;
        }

#if PVP_KINETIC_ONLY
      tmp =  gshell1->coefficient_unnorm(gc1,i)
             * gshell2->coefficient_unnorm(gc2,j);
#else
      // The "-" comes from i*i in the momentum operators
      tmp =  -gshell1->coefficient_unnorm(gc1,i)
             * gshell2->coefficient_unnorm(gc2,j)
             * coef;
#endif

      double *primbuffer;
      int i2,j2,k2;
      int tsize;

      //////////////////////////////////////////
      // Compute phi_1 V del^2 phi_2

      // differentiate the cartesian prefactors twice:
      if (b >= 2) {
#if PVP_KINETIC_ONLY
        primbuffer = inter(a,b-2,0);
        ss =   pow(M_PI/(gshell1->exponent(i)
                         +gshell2->exponent(j)),1.5)
               * exp(- oozeta * gshell1->exponent(i)
                     * gshell2->exponent(j) * AmB2);
        int i1,j1,k1;
        FOR_CART(i1,j1,k1,a)
          FOR_CART(i2,j2,k2,b-2)
            *primbuffer++ = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
          END_FOR_CART
        END_FOR_CART
#else
        comp_prim_block_nuclear(a,b-2);
#endif

        primbuffer = inter(a,b-2,0);
        tsize = INT_NCART(b-2);
        for (int ip=0; ip<size1; ip++) {
          int jp=0;
          FOR_CART(i2,j2,k2,b)
            if (i2 >= 2) {
              int tjp = INT_CARTINDEX(b-2,i2-2,j2);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                += i2 * (i2-1) * tmp * primbuffer[tjp];
#if PVP_DEBUG
              std::cout << "contribA1x = "
                        << i2 * (i2-1) * tmp * primbuffer[tjp]
                        << std::endl;
#endif
              }
            if (j2 >= 2) {
              int tjp = INT_CARTINDEX(b-2,i2,j2-2);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                += j2 * (j2-1) * tmp * primbuffer[tjp];
#if PVP_DEBUG
              std::cout << "contribA1y = "
                        << j2 * (j2-1) * tmp * primbuffer[tjp]
                        << std::endl;
#endif
              }
            if (k2 >= 2) {
              int tjp = INT_CARTINDEX(b-2,i2,j2);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                += k2 * (k2-1) * tmp * primbuffer[tjp];
#if PVP_DEBUG
              std::cout << "contribA1z = "
                        << k2 * (k2-1) * tmp * primbuffer[tjp]
                        << std::endl;
#endif
              }
            jp++;
          END_FOR_CART
          primbuffer += tsize;
        }

      }

      // differentiate the cartesian prefactor and the exponential
#if PVP_KINETIC_ONLY
      primbuffer = inter(a,b,0);
      ss =   pow(M_PI/(gshell1->exponent(i)
                       +gshell2->exponent(j)),1.5)
             * exp(- oozeta * gshell1->exponent(i)
                   * gshell2->exponent(j) * AmB2);
      int i1,j1,k1;
      FOR_CART(i1,j1,k1,a)
        FOR_CART(i2,j2,k2,b)
          *primbuffer++ = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
        END_FOR_CART
      END_FOR_CART
#if PVP_DEBUG
      std::cout << "primbuffer[0] (overlap a b) = " << *inter(a,b,0)
                << ", " << "*tmp = " << *inter(a,b,0) * tmp
                << std::endl;
#endif
#else
      comp_prim_block_nuclear(a,b);
#endif

      primbuffer = inter(a,b,0);
      for (int ip=0; ip<size1; ip++) {
        int jp=0;
        FOR_CART(i2,j2,k2,b)
          double val = *primbuffer;
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            -= 2.0 * (2.0*i2+1.0) * alpha2 * tmp * val;
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            -= 2.0 * (2.0*j2+1.0) * alpha2 * tmp * val;
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            -= 2.0 * (2.0*k2+1.0) * alpha2 * tmp * val;
#if PVP_DEBUG
          std::cout << "contribA2x = "
                    << - 2.0 * (2.0*i2+1.0) * alpha2 * tmp * val
                    << std::endl;
          std::cout << "contribA2y = "
                    << - 2.0 * (2.0*j2+1.0) * alpha2 * tmp * val
                    << std::endl;
          std::cout << "contribA2z = "
                    << - 2.0 * (2.0*k2+1.0) * alpha2 * tmp * val
                    << std::endl;
#endif
          jp++;
          primbuffer++;
        END_FOR_CART
      }

      // differentiate exponential twice
#if PVP_KINETIC_ONLY
      primbuffer = inter(a,b+2,0);
      ss =   pow(M_PI/(gshell1->exponent(i)
                       +gshell2->exponent(j)),1.5)
             * exp(- oozeta * gshell1->exponent(i)
                   * gshell2->exponent(j) * AmB2);
      FOR_CART(i1,j1,k1,a)
        FOR_CART(i2,j2,k2,b+2)
          *primbuffer++ = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
        END_FOR_CART
      END_FOR_CART
#if PVP_DEBUG
      std::cout << "primbuffer[0] (overlap a b+2) = " << *inter(a,b+2,0)
                << std::endl;
#endif
#else
      comp_prim_block_nuclear(a,b+2);
#endif

      primbuffer = inter(a,b+2,0);
      tsize = INT_NCART(b+2);
      for (int ip=0; ip<size1; ip++) {
        int jp=0;
        FOR_CART(i2,j2,k2,b)
          int tjp_x = INT_CARTINDEX(b+2,i2+2,j2);
          int tjp_y = INT_CARTINDEX(b+2,i2,j2+2);
          int tjp_z = INT_CARTINDEX(b+2,i2,j2);
          double val = primbuffer[tjp_x]
                       +primbuffer[tjp_y]
                       +primbuffer[tjp_z];
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            += 4.0 * alpha2 * alpha2 * tmp * val;
#if PVP_DEBUG
          std::cout << "contribA3 = "
                    << 4.0 * alpha2 * alpha2 * tmp * val
                    << std::endl;
#endif
          jp++;
        END_FOR_CART
        primbuffer += tsize;
      }

      //////////////////////////////////////////
      // Compute phi_1 (del V) del phi_2

#if !PVP_KINETIC_ONLY

      // differentiate the cartesian prefactor
      if (b >= 1) {
        comp_prim_block_efield(a,b-1);
        primbuffer = efield_inter(a,b-1,0);
        int tsize = INT_NCART(b-1);
        for (int ip=0; ip<size1; ip++) {
          int jp=0;
          FOR_CART(i2,j2,k2,b)
            if (i2 >= 1) {
              int tjp = INT_CARTINDEX(b-1,i2-1,j2);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                -= i2 * tmp * primbuffer[3*tjp+0];
#if PVP_DEBUG
              std::cout << "contribB1x = "
                        << - i2 * tmp * primbuffer[3*tjp+0]
                        << endl;
#endif
              }
            if (j2 >= 1) {
              int tjp = INT_CARTINDEX(b-1,i2,j2-1);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                -= j2 * tmp * primbuffer[3*tjp+1];
#if PVP_DEBUG
              std::cout << "contribB1y = "
                        << - j2 * tmp * primbuffer[3*tjp+1]
                        << endl;
#endif
              }
            if (k2 >= 1) {
              int tjp = INT_CARTINDEX(b-1,i2,j2);
              buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
                -= k2 * tmp * primbuffer[3*tjp+2];
#if PVP_DEBUG
              std::cout << "contribB1z = "
                        << - k2 * tmp * primbuffer[3*tjp+2]
                        << endl;
#endif
              }
            jp++;
          END_FOR_CART
          primbuffer += 3*tsize;
        }
      }

      // differentiate the exponential
      comp_prim_block_efield(a,b+1);
      primbuffer = efield_inter(a,b+1,0);
      tsize = INT_NCART(b+1);
      for (int ip=0; ip<size1; ip++) {
        int jp=0;
        FOR_CART(i2,j2,k2,b)
          int tip = INT_CARTINDEX(b+1,i2+1,j2);
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            += 2.0 * alpha2 * tmp * primbuffer[3*tip+0];

          int tjp = INT_CARTINDEX(b+1,i2,j2+1);
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            += 2.0 * alpha2 * tmp * primbuffer[3*tjp+1];

          int tkp = INT_CARTINDEX(b+1,i2,j2);
          buffer[(ip+gcoff1)*gcsize2+jp+gcoff2]
            += 2.0 * alpha2 * tmp * primbuffer[3*tkp+2];

#if PVP_DEBUG
          std::cout << "contribB2x = "
                    << 2.0 * alpha2 * tmp * primbuffer[3*tip+0]
                    << endl;

          std::cout << "contribB2y = "
                    << 2.0 * alpha2 * tmp * primbuffer[3*tjp+1]
                    << endl;

          std::cout << "contribB2z = "
                    << 2.0 * alpha2 * tmp * primbuffer[3*tkp+2]
                    << endl;
#endif

          jp++;
        END_FOR_CART
        primbuffer += 3*tsize;
        }

#endif


      }
    }
  }

void
Int1eV3::comp_prim_block_nuclear(int a, int b)
{
  int im, ia, ib;
  int l = a + b;

  // fill in the ia+ib=0 integrals
  for (im=0; im<=l; im++) {
#if DEBUG_NUC_PRIM > 1
    ExEnv::outn() << "BUILD: M = " << im
         << " A = " << 0
         << " B = " << 0
         << endl;
#endif
    inter(0,0,im)[0] = fjttable_[im];
    }

  for (im=l-1; im>=0; im--) {
    int lm = l-im;
    // build the integrals for a = 0
    for (ib=1; ib<=lm && ib<=b; ib++) {
#if DEBUG_NUC_PRIM > 1
      ExEnv::outn() << "BUILD: M = " << im
           << " A = " << 0
           << " B = " << ib
           << endl;
#endif
      comp_prim_block_nuclear_build_b(ib,im);
      }
    for (ia=1; ia<=lm && ia<=a; ia++) {
      for (ib=0; ib<=lm-ia && ib<=b; ib++) {
#if DEBUG_NUC_PRIM > 1
        ExEnv::outn() << "BUILD: M = " << im
             << " A = " << ia
             << " B = " << ib
             << endl;
#endif
        comp_prim_block_nuclear_build_a(ia,ib,im);
        }
      }
    }

#if DEBUG_NUC_PRIM
  for (im=0; im<=l; im++) {
    int lm = l-im;
    for (ia=0; ia<=lm && ia<=a; ia++) {
      int na = INT_NCART_NN(a);
      for (ib=0; ib<=lm-ia && ib<=b; ib++) {
        int nb = INT_NCART_NN(b);
#if DEBUG_NUC_PRIM > 1
        ExEnv::outn() << "M = " << im
             << " A = " << ia
             << " B = " << ib
             << endl;
#endif
        double *buf = inter(ia,ib,im);
        int i1,j1,k1, i2,j2,k2;
        FOR_CART(i1,j1,k1,ia) {
          FOR_CART(i2,j2,k2,ib) {
            double fast = *buf++;
            double slow = comp_prim_nuclear(i1, j1, k1,
                                            i2, j2, k2, im);
            if (fast > 999.0) fast = 999.0;
            if (fast < -999.0) fast = -999.0;
            if (fabs(fast-slow)>1.0e-12) {
              ExEnv::outn() << scprintf("(%d %d %d) (%d %d %d) (%d): ",
                               i1,j1,k1, i2,j2,k2, im)
                   << scprintf(" % 20.15f % 20.15f",
                               fast, slow)
                   << endl;
              }
            } END_FOR_CART;
          } END_FOR_CART;
        }
      }
    }
#endif
}

void
Int1eV3::comp_prim_block_nuclear_build_a(int a, int b, int m)
{
  double *I000 = inter(a,b,m);
  double *I100 = inter(a-1,b,m);
  double *I101 = inter(a-1,b,m+1);
  double *I200 = (a>1?inter(a-2,b,m):0);
  double *I201 = (a>1?inter(a-2,b,m+1):0);
  double *I110 = (b?inter(a-1,b-1,m):0);
  double *I111 = (b?inter(a-1,b-1,m+1):0);
  int i1,j1,k1;
  int i2,j2,k2;
  int carta=0;
  int sizeb = INT_NCART_NN(b);
  int sizebm1 = INT_NCART_DEC(b,sizeb);
  FOR_CART(i1,j1,k1,a) {
    int cartb=0;
    FOR_CART(i2,j2,k2,b) {
      double result = 0.0;
      if (i1) {
        int am1 = INT_CARTINDEX(a-1,i1-1,j1);
        result  = PmA[0] * I100[am1*sizeb+cartb];
        result -= PmC[0] * I101[am1*sizeb+cartb];
        if (i1>1) {
          int am2 = INT_CARTINDEX(a-2,i1-2,j1);
          result += oo2zeta * (i1-1)
                    *(I200[am2*sizeb+cartb] - I201[am2*sizeb+cartb]);
          }
        if (i2) {
          int bm1 = INT_CARTINDEX(b-1,i2-1,j2);
          result += oo2zeta * i2
                    *(I110[am1*sizebm1+bm1] - I111[am1*sizebm1+bm1]);
          }
        }
      else if (j1) {
        int am1 = INT_CARTINDEX(a-1,i1,j1-1);
        result  = PmA[1] * I100[am1*sizeb+cartb];
        result -= PmC[1] * I101[am1*sizeb+cartb];
        if (j1>1) {
          int am2 = INT_CARTINDEX(a-2,i1,j1-2);
          result += oo2zeta * (j1-1)
                    *(I200[am2*sizeb+cartb] - I201[am2*sizeb+cartb]);
          }
        if (j2) {
          int bm1 = INT_CARTINDEX(b-1,i2,j2-1);
          result += oo2zeta * j2
                    *(I110[am1*sizebm1+bm1] - I111[am1*sizebm1+bm1]);
          }
        }
      else if (k1) {
        int am1 = INT_CARTINDEX(a-1,i1,j1);
        result  = PmA[2] * I100[am1*sizeb+cartb];
        result -= PmC[2] * I101[am1*sizeb+cartb];
        if (k1>1) {
          int am2 = INT_CARTINDEX(a-2,i1,j1);
          result += oo2zeta * (k1-1)
                    *(I200[am2*sizeb+cartb] - I201[am2*sizeb+cartb]);
          }
        if (k2) {
          int bm1 = INT_CARTINDEX(b-1,i2,j2);
          result += oo2zeta * k2
                    *(I110[am1*sizebm1+bm1] - I111[am1*sizebm1+bm1]);
          }
        }
      I000[carta*sizeb+cartb] = result;
      cartb++;
      } END_FOR_CART;
    carta++;
    } END_FOR_CART;
}

void
Int1eV3::comp_prim_block_nuclear_build_b(int b, int m)
{
  double *I000 = inter(0,b,m);
  double *I010 = inter(0,b-1,m);
  double *I011 = inter(0,b-1,m+1);
  double *I020 = (b>1?inter(0,b-2,m):0);
  double *I021 = (b>1?inter(0,b-2,m+1):0);
  int i2,j2,k2;
  int cartb=0;
  FOR_CART(i2,j2,k2,b) {
    double result = 0.0;

    if (i2) {
      int bm1 = INT_CARTINDEX(b-1,i2-1,j2);
      result  = PmB[0] * I010[bm1];
      result -= PmC[0] * I011[bm1];
      if (i2>1) {
        int bm2 = INT_CARTINDEX(b-2,i2-2,j2);
        result += oo2zeta * (i2-1) * (I020[bm2] - I021[bm2]);
        }
      }
    else if (j2) {
      int bm1 = INT_CARTINDEX(b-1,i2,j2-1);
      result  = PmB[1] * I010[bm1];
      result -= PmC[1] * I011[bm1];
      if (j2>1) {
        int bm2 = INT_CARTINDEX(b-2,i2,j2-2);
        result += oo2zeta * (j2-1) * (I020[bm2] - I021[bm2]);
        }
      }
    else if (k2) {
      int bm1 = INT_CARTINDEX(b-1,i2,j2);
      result  = PmB[2] * I010[bm1];
      result -= PmC[2] * I011[bm1];
      if (k2>1) {
        int bm2 = INT_CARTINDEX(b-2,i2,j2);
        result += oo2zeta * (k2-1) * (I020[bm2] - I021[bm2]);
        }
      }

    I000[cartb] = result;
    cartb++;
    } END_FOR_CART;
}

double
Int1eV3::comp_prim_nuclear(int i1, int j1, int k1,
                           int i2, int j2, int k2, int m)
{
  double result;

  if (i1) {
    result  = PmA[0] * comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m);
    result -= PmC[0] * comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m+1);
    if (i1>1) result += oo2zeta * (i1-1)
                       * (  comp_prim_nuclear(i1-2,j1,k1,i2,j2,k2,m)
                          - comp_prim_nuclear(i1-2,j1,k1,i2,j2,k2,m+1));
    if (i2) result += oo2zeta * i2
                     * (  comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m+1));
    }
  else if (j1) {
    result  = PmA[1] * comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m);
    result -= PmC[1] * comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m+1);
    if (j1>1) result += oo2zeta * (j1-1)
                       * (  comp_prim_nuclear(i1,j1-2,k1,i2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1-2,k1,i2,j2,k2,m+1));
    if (j2) result += oo2zeta * j2
                     * (  comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m+1));
    }
  else if (k1) {
    result  = PmA[2] * comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m);
    result -= PmC[2] * comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m+1);
    if (k1>1) result += oo2zeta * (k1-1)
                       * (  comp_prim_nuclear(i1,j1,k1-2,i2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1-2,i2,j2,k2,m+1));
    if (k2) result += oo2zeta * k2
                     * (  comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m+1));
    }
  else if (i2) {
    result  = PmB[0] * comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m);
    result -= PmC[0] * comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m+1);
    if (i2>1) result += oo2zeta * (i2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2-2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2-2,j2,k2,m+1));
    if (i1) result += oo2zeta * i1
                     * (  comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m+1));
    }
  else if (j2) {
    result  = PmB[1] * comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m);
    result -= PmC[1] * comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m+1);
    if (j2>1) result += oo2zeta * (j2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2,j2-2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2,j2-2,k2,m+1));
    if (j1) result += oo2zeta * j1
                     * (  comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m+1));
    }
  else if (k2) {
    result  = PmB[2] * comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m);
    result -= PmC[2] * comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m+1);
    if (k2>1) result += oo2zeta * (k2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2,j2,k2-2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2,j2,k2-2,m+1));
    if (k1) result += oo2zeta * k1
                     * (  comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m+1));
    }
  else result = fjttable_[m];

  return result;
  }

/* Compute the electric field integral for the shell.  The arguments are the
 * the electric field vector, the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
void
Int1eV3::comp_shell_efield(double *efield,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2)
{
  int i,j,k,xyz;
  double result[3];
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  int am;

  am = i1+j1+k1+i2+j2+k2;

  /* Loop over the primitives in the shells. */
  for (xyz=0; xyz<3; xyz++) result[xyz] = 0.0;
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      zeta = gshell1->exponent(i) + gshell2->exponent(j);
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz] + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * M_PI/(gshell1->exponent(i)
                                           +gshell2->exponent(j))
           * exp(- oozeta * gshell1->exponent(i)
                 * gshell2->exponent(j) * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am+1,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am+1; k++) {
        fjttable_[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      for (xyz=0; xyz<3; xyz++) {
        result[xyz] +=  gshell1->coefficient_unnorm(gc1,i)
                      * gshell2->coefficient_unnorm(gc2,j)
                      * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2,0);
        }
      }
    }

  for (xyz=0; xyz<3; xyz++) efield[xyz] = result[xyz];

  }

/* Compute the electric field integral for the shell.  The arguments are the
 * the electric field vector, the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
void
Int1eV3::comp_shell_block_efield(int gc1, int a, int gc2, int b,
                                 int gcsize2, int gcoff1, int gcoff2,
                                 double coef, double *buffer)
{
  int i,j,k,xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  int am = a + b;
  int size1 = INT_NCART_NN(a);
  int size2 = INT_NCART_NN(b);

  /* Loop over the primitives in the shells. */
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      zeta = gshell1->exponent(i) + gshell2->exponent(j);
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(gshell1->exponent(i) * A[xyz]
                     + gshell2->exponent(j) * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * M_PI/(gshell1->exponent(i)
                                           +gshell2->exponent(j))
           * exp(- oozeta * gshell1->exponent(i)
                 * gshell2->exponent(j) * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am+1,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am+1; k++) {
        fjttable_[k] *= auxcoef;
        }

      double tmp = gshell1->coefficient_unnorm(gc1,i)
                   * gshell2->coefficient_unnorm(gc2,j)
                   * coef;

      /* Compute the nuclear attraction integral. */
      comp_prim_block_efield(a,b);

      double *primbuffer = efield_inter(a,b,0);
      for (int ip=0; ip<size1; ip++) {
        for (int jp=0; jp<size2; jp++) {
          for (xyz=0; xyz<3; xyz++) {
            buffer[((ip+gcoff1)*gcsize2+jp+gcoff2)*3+xyz]
              += tmp * *primbuffer++;
            }
          }
        }
      }
    }

  }

void
Int1eV3::comp_prim_block_efield(int a, int b)
{
  int  xyz;
  int im, ia, ib;
  int l = a + b;

  // Fill in the nuclear integrals which are needed as intermediates.
  // m=0 is not needed.

  // fill in the ia+ib=0 integrals, skipping m=0
  for (im=1; im<=l; im++) {
#if DEBUG_EFIELD_PRIM > 1
    ExEnv::outn() << "BUILD NUC: M = " << im
         << " A = " << 0
         << " B = " << 0
         << endl;
#endif
    inter(0,0,im)[0] = fjttable_[im];
    }

  // skipping m=0
  for (im=l-1; im>0; im--) {
    int lm = l-im;
    // build the integrals for a = 0
    for (ib=1; ib<=lm && ib<=b; ib++) {
#if DEBUG_EFIELD_PRIM > 1
      ExEnv::outn() << "BUILD NUC: M = " << im
           << " A = " << 0
           << " B = " << ib
           << endl;
#endif
      comp_prim_block_nuclear_build_b(ib,im);
      }
    for (ia=1; ia<=lm && ia<=a; ia++) {
      for (ib=0; ib<=lm-ia && ib<=b; ib++) {
#if DEBUG_EFIELD_PRIM > 1
        ExEnv::outn() << "BUILD NUC: M = " << im
             << " A = " << ia
             << " B = " << ib
             << endl;
#endif
        comp_prim_block_nuclear_build_a(ia,ib,im);
        }
      }
    }

  // now complete the efield integrals

  // fill in the ia+ib=0 integrals
  for (im=0; im<=l; im++) {
#if DEBUG_EFIELD_PRIM > 1
    ExEnv::outn() << "BUILD EFIELD: M = " << im
         << " A = " << 0
         << " B = " << 0
         << endl;
#endif
    double *tmp = efield_inter(0,0,im);
    for (xyz=0; xyz<3; xyz++) {
      *tmp++ = 2.0 * zeta * PmC[xyz] * fjttable_[im+1];
      }
    }

  for (im=l-1; im>=0; im--) {
    int lm = l-im;
    // build the integrals for a = 0
    for (ib=1; ib<=lm && ib<=b; ib++) {
#if DEBUG_EFIELD_PRIM > 1
      ExEnv::outn() << "BUILD EFIELD: M = " << im
           << " A = " << 0
           << " B = " << ib
           << endl;
#endif
      comp_prim_block_efield_build_b(ib,im);
      }
    for (ia=1; ia<=lm && ia<=a; ia++) {
      for (ib=0; ib<=lm-ia && ib<=b; ib++) {
#if DEBUG_EFIELD_PRIM > 1
        ExEnv::outn() << "BUILD EFIELD: M = " << im
             << " A = " << ia
             << " B = " << ib
             << endl;
#endif
        comp_prim_block_efield_build_a(ia,ib,im);
        }
      }
    }

#if DEBUG_EFIELD_PRIM
  for (im=0; im<=l; im++) {
    int lm = l-im;
    for (ia=0; ia<=lm && ia<=a; ia++) {
      int na = INT_NCART_NN(a);
      for (ib=0; ib<=lm-ia && ib<=b; ib++) {
        int nb = INT_NCART_NN(b);
#if DEBUG_EFIELD_PRIM > 1
        ExEnv::outn() << "M = " << im
             << " A = " << ia
             << " B = " << ib
             << endl;
#endif
        double *buf = efield_inter(ia,ib,im);
        int i1,j1,k1, i2,j2,k2;
        FOR_CART(i1,j1,k1,ia) {
          FOR_CART(i2,j2,k2,ib) {
            for (xyz=0; xyz<3; xyz++) {
              double fast = *buf++;
              double slow = comp_prim_efield(xyz, i1, j1, k1,
                                             i2, j2, k2, im);
              if (fast > 999.0) fast = 999.0;
              if (fast < -999.0) fast = -999.0;
              if (fabs(fast-slow)>1.0e-12) {
                ExEnv::outn() << scprintf("(%d %d %d) (%d %d %d) (%d): ",
                                 i1,j1,k1, i2,j2,k2, im)
                     << scprintf(" % 20.15f % 20.15f",
                                 fast, slow)
                     << endl;
                }
              }
            } END_FOR_CART;
          } END_FOR_CART;
        }
      }
    }
#endif
}

void
Int1eV3::comp_prim_block_efield_build_a(int a, int b, int m)
{
  double *I000 = efield_inter(a,b,m);
  double *I100 = efield_inter(a-1,b,m);
  double *I101 = efield_inter(a-1,b,m+1);
  double *I200 = (a>1?efield_inter(a-2,b,m):0);
  double *I201 = (a>1?efield_inter(a-2,b,m+1):0);
  double *I110 = (b?efield_inter(a-1,b-1,m):0);
  double *I111 = (b?efield_inter(a-1,b-1,m+1):0);
  double *nucI101 = inter(a-1,b,m+1);
  int i1,j1,k1;
  int i2,j2,k2;
  int xyz;
  int carta=0;
  int sizeb = INT_NCART_NN(b);
  int sizebm1 = INT_NCART_DEC(b,sizeb);
  FOR_CART(i1,j1,k1,a) {
    int cartb=0;
    FOR_CART(i2,j2,k2,b) {
      for (xyz=0; xyz<3; xyz++) {
        double result = 0.0;
        if (i1) {
          int am1 = INT_CARTINDEX(a-1,i1-1,j1);
          result  = PmA[0] * I100[(am1*sizeb+cartb)*3+xyz];
          result -= PmC[0] * I101[(am1*sizeb+cartb)*3+xyz];
          if (i1>1) {
            int am2 = INT_CARTINDEX(a-2,i1-2,j1);
            result += oo2zeta * (i1-1)
                      *(I200[(am2*sizeb+cartb)*3+xyz]
                        - I201[(am2*sizeb+cartb)*3+xyz]);
            }
          if (i2) {
            int bm1 = INT_CARTINDEX(b-1,i2-1,j2);
            result += oo2zeta * i2
                      *(I110[(am1*sizebm1+bm1)*3+xyz]
                        - I111[(am1*sizebm1+bm1)*3+xyz]);
            }
          if (xyz==0) result += nucI101[am1*sizeb+cartb];
          }
        else if (j1) {
          int am1 = INT_CARTINDEX(a-1,i1,j1-1);
          result  = PmA[1] * I100[(am1*sizeb+cartb)*3+xyz];
          result -= PmC[1] * I101[(am1*sizeb+cartb)*3+xyz];
          if (j1>1) {
            int am2 = INT_CARTINDEX(a-2,i1,j1-2);
            result += oo2zeta * (j1-1)
                      *(I200[(am2*sizeb+cartb)*3+xyz]
                        - I201[(am2*sizeb+cartb)*3+xyz]);
            }
          if (j2) {
            int bm1 = INT_CARTINDEX(b-1,i2,j2-1);
            result += oo2zeta * j2
                      *(I110[(am1*sizebm1+bm1)*3+xyz]
                        - I111[(am1*sizebm1+bm1)*3+xyz]);
            }
          if (xyz==1) result += nucI101[am1*sizeb+cartb];
          }
        else if (k1) {
          int am1 = INT_CARTINDEX(a-1,i1,j1);
          result  = PmA[2] * I100[(am1*sizeb+cartb)*3+xyz];
          result -= PmC[2] * I101[(am1*sizeb+cartb)*3+xyz];
          if (k1>1) {
            int am2 = INT_CARTINDEX(a-2,i1,j1);
            result += oo2zeta * (k1-1)
                      *(I200[(am2*sizeb+cartb)*3+xyz]
                        - I201[(am2*sizeb+cartb)*3+xyz]);
            }
          if (k2) {
            int bm1 = INT_CARTINDEX(b-1,i2,j2);
            result += oo2zeta * k2
                      *(I110[(am1*sizebm1+bm1)*3+xyz]
                        - I111[(am1*sizebm1+bm1)*3+xyz]);
            }
          if (xyz==2) result += nucI101[am1*sizeb+cartb];
          }
        I000[(carta*sizeb+cartb)*3+xyz] = result;
        }
      cartb++;
      } END_FOR_CART;
    carta++;
    } END_FOR_CART;
}

void
Int1eV3::comp_prim_block_efield_build_b(int b, int m)
{
  double *I000 = efield_inter(0,b,m);
  double *I010 = efield_inter(0,b-1,m);
  double *I011 = efield_inter(0,b-1,m+1);
  double *I020 = (b>1?efield_inter(0,b-2,m):0);
  double *I021 = (b>1?efield_inter(0,b-2,m+1):0);
  double *nucI011 = inter(0,b-1,m+1);
  int xyz;
  int i2,j2,k2;
  int cartb=0;
  FOR_CART(i2,j2,k2,b) {
    for (xyz=0; xyz<3; xyz++) {
      double result = 0.0;

      if (i2) {
        int bm1 = INT_CARTINDEX(b-1,i2-1,j2);
        result  = PmB[0] * I010[(bm1)*3+xyz];
        result -= PmC[0] * I011[(bm1)*3+xyz];
        if (i2>1) {
          int bm2 = INT_CARTINDEX(b-2,i2-2,j2);
          result += oo2zeta * (i2-1) * (I020[(bm2)*3+xyz]
                                        - I021[(bm2)*3+xyz]);
          }
        if (xyz==0) result += nucI011[bm1];
        }
      else if (j2) {
        int bm1 = INT_CARTINDEX(b-1,i2,j2-1);
        result  = PmB[1] * I010[(bm1)*3+xyz];
        result -= PmC[1] * I011[(bm1)*3+xyz];
        if (j2>1) {
          int bm2 = INT_CARTINDEX(b-2,i2,j2-2);
          result += oo2zeta * (j2-1) * (I020[(bm2)*3+xyz]
                                        - I021[(bm2)*3+xyz]);
          }
        if (xyz==1) result += nucI011[bm1];
        }
      else if (k2) {
        int bm1 = INT_CARTINDEX(b-1,i2,j2);
        result  = PmB[2] * I010[(bm1)*3+xyz];
        result -= PmC[2] * I011[(bm1)*3+xyz];
        if (k2>1) {
          int bm2 = INT_CARTINDEX(b-2,i2,j2);
          result += oo2zeta * (k2-1) * (I020[(bm2)*3+xyz]
                                        - I021[(bm2)*3+xyz]);
          }
        if (xyz==2) result += nucI011[bm1];
        }

      I000[(cartb)*3+xyz] = result;
      }
    cartb++;
    } END_FOR_CART;
}

double
Int1eV3::comp_prim_efield(int xyz, int i1, int j1, int k1,
                          int i2, int j2, int k2, int m)
{
  double result;

  /* if ((xyz != 0) || (i1 != 1)) return 0.0; */

  if (i1) {
    result  = PmA[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m);
    result -= PmC[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m+1);
    if (i1>1) result += oo2zeta * (i1-1)
                       * (  comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m+1));
    if (i2) result += oo2zeta * i2
                     * (  comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m+1));
    if (xyz==0) result += comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m+1);
    }
  else if (j1) {
    result  = PmA[1] * comp_prim_efield(xyz,i1,j1-1,k1,i2,j2,k2,m);
    result -= PmC[1] * comp_prim_efield(xyz,i1,j1-1,k1,i2,j2,k2,m+1);
    if (j1>1) result += oo2zeta * (j1-1)
                       * (  comp_prim_efield(xyz,i1,j1-2,k1,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1-2,k1,i2,j2,k2,m+1));
    if (j2) result += oo2zeta * j2
                     * (  comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m+1));
    if (xyz==1) result += comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m+1);
    }
  else if (k1) {
    result  = PmA[2] * comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2,m);
    result -= PmC[2] * comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2,m+1);
    if (k1>1) result += oo2zeta * (k1-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1-2,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1-2,i2,j2,k2,m+1));
    if (k2) result += oo2zeta * k2
                     * (  comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m+1));
    if (xyz==2) result += comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m+1);
    }
  else if (i2) {
    result  = PmB[0] * comp_prim_efield(xyz,i1,j1,k1,i2-1,j2,k2,m);
    result -= PmC[0] * comp_prim_efield(xyz,i1,j1,k1,i2-1,j2,k2,m+1);
    if (i2>1) result += oo2zeta * (i2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2-2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2-2,j2,k2,m+1));
    if (i1) result += oo2zeta * i1
                     * (  comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m+1));
    if (xyz==0) result += comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m+1);
    }
  else if (j2) {
    result  = PmB[1] * comp_prim_efield(xyz,i1,j1,k1,i2,j2-1,k2,m);
    result -= PmC[1] * comp_prim_efield(xyz,i1,j1,k1,i2,j2-1,k2,m+1);
    if (j2>1) result += oo2zeta * (j2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2,j2-2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2,j2-2,k2,m+1));
    if (j1) result += oo2zeta * j1
                     * (  comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m+1));
    if (xyz==1) result += comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m+1);
    }
  else if (k2) {
    result  = PmB[2] * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-1,m);
    result -= PmC[2] * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-1,m+1);
    if (k2>1) result += oo2zeta * (k2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-2,m+1));
    if (k1) result += oo2zeta * k1
                     * (  comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m+1));
    if (xyz==2) result += comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m+1);
    }
  else {
    /* We arrive here if we have a (s| |s) type efield integral.
     * The fjttable contains the standard (s| |s) nuc attr integrals.
     */
    result = 2.0 * zeta * PmC[xyz] * fjttable_[m+1];
    }

  return result;
  }


/* --------------------------------------------------------------- */
/* ------------- Routines for dipole moment integrals ------------ */
/* --------------------------------------------------------------- */

/* This computes the dipole integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The last arg, com, is the origin of the coordinate
 * system used to compute the dipole moment.
 */
void
Int1eV3::int_accum_shell_dipole(int ish, int jsh,
                                double *com)
{
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int gc1,gc2;
  int index,index1,index2;
  double dipole[3];
  int xyz;

  for (xyz=0; xyz<3; xyz++) C[xyz] = com[xyz];

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;
  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2);
      for(mu=0; mu < 3; mu++) {
        cartesianbuffer[index] = dipole[mu];
        index++;
        }
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)
  accum_transform_1e_xyz(integral_,
                         cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the dipole integrals between functions in two shells.
 * The result is placed in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The last arg, com, is the origin of the coordinate
 * system used to compute the dipole moment.
 */
void
Int1eV3::dipole(int ish, int jsh, double *com)
{
  int c1,i1,j1,k1,c2,i2,j2,k2;
  int gc1,gc2;
  int index,index1,index2;
  double dipole[3];
  int xyz;

  for (xyz=0; xyz<3; xyz++) C[xyz] = com[xyz];

  c1 = bs1_->shell_to_center(ish);
  c2 = bs2_->shell_to_center(jsh);
  for (xyz=0; xyz<3; xyz++) {
      A[xyz] = bs1_->r(c1,xyz);
      B[xyz] = bs2_->r(c2,xyz);
    }
  gshell1 = &bs1_->shell(ish);
  gshell2 = &bs2_->shell(jsh);
  index = 0;
  FOR_GCCART_GS(gc1,index1,i1,j1,k1,gshell1)
    FOR_GCCART_GS(gc2,index2,i2,j2,k2,gshell2)
      comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2);
      for(mu=0; mu < 3; mu++) {
        cartesianbuffer[index] = dipole[mu];
        index++;
        }
      END_FOR_GCCART_GS(index2)
    END_FOR_GCCART_GS(index1)
  transform_1e_xyz(integral_,
                   cartesianbuffer, buff, gshell1, gshell2);
  }

void
Int1eV3::comp_shell_dipole(double* dipole,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2)
{
  double exp1,exp2;
  int i,j,xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  dipole[0] = dipole[1] = dipole[2] = 0.0;

  if ((i1<0)||(j1<0)||(k1<0)||(i2<0)||(j2<0)||(k2<0)) return;

  /* Loop over the primitives in the shells. */
  for (i=0; i<gshell1->nprimitive(); i++) {
    for (j=0; j<gshell2->nprimitive(); j++) {

      /* Compute the intermediates. */
      exp1 = gshell1->exponent(i);
      exp2 = gshell2->exponent(j);
      oozeta = 1.0/(exp1 + exp2);
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(exp1 * A[xyz] + exp2 * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      ss =   pow(M_PI/(exp1+exp2),1.5)
           * exp(- oozeta * exp1 * exp2 * AmB2);
      for (mu=0; mu<3; mu++) sMus[mu] = ss * PmC[mu];
      tmp     =  gshell1->coefficient_unnorm(gc1,i)
               * gshell2->coefficient_unnorm(gc2,j);
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      dipole[0] += tmp * comp_prim_dipole(0,i1,j1,k1,i2,j2,k2);
      dipole[1] += tmp * comp_prim_dipole(1,i1,j1,k1,i2,j2,k2);
      dipole[2] += tmp * comp_prim_dipole(2,i1,j1,k1,i2,j2,k2);
      }
    }

  }

double
Int1eV3::comp_prim_dipole(int axis,
                          int i1, int j1, int k1,
                          int i2, int j2, int k2)
{
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_dipole(axis,i1-1,j1,k1,i2,j2,k2);
    if (i2) 
      result += oo2zeta*i2*comp_prim_dipole(axis,i1-1,j1,k1,i2-1,j2,k2);
    if (i1>1)
      result += oo2zeta*(i1-1)*comp_prim_dipole(axis,i1-2,j1,k1,i2,j2,k2);
    if(axis==0) result += oo2zeta*comp_prim_overlap(i1-1,j1,k1,i2,j2,k2);
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_dipole(axis,i1,j1-1,k1,i2,j2,k2);
    if (j2) 
      result += oo2zeta*j2*comp_prim_dipole(axis,i1,j1-1,k1,i2,j2-1,k2);
    if (j1>1) 
      result += oo2zeta*(j1-1)*comp_prim_dipole(axis,i1,j1-2,k1,i2,j2,k2);
    if(axis==1) result += oo2zeta*comp_prim_overlap(i1,j1-1,k1,i2,j2,k2);
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_dipole(axis,i1,j1,k1-1,i2,j2,k2);
    if (k2) 
      result += oo2zeta*k2*comp_prim_dipole(axis,i1,j1,k1-1,i2,j2,k2-1);
    if (k1>1) 
      result += oo2zeta*(k1-1)*comp_prim_dipole(axis,i1,j1,k1-2,i2,j2,k2);
    if(axis==2) result += oo2zeta*comp_prim_overlap(i1,j1,k1-1,i2,j2,k2);
    return result;
    }
  if (i2) {
    result = PmB[0] * comp_prim_dipole(axis,i1,j1,k1,i2-1,j2,k2);
    if (i1) 
      result += oo2zeta*i1*comp_prim_dipole(axis,i1-1,j1,k1,i2-1,j2,k2);
    if (i2>1) 
      result += oo2zeta*(i2-1)*comp_prim_dipole(axis,i1,j1,k1,i2-2,j2,k2);
    if(axis==0) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_dipole(axis,i1,j1,k1,i2,j2-1,k2);
    if (j1) 
      result += oo2zeta*i1*comp_prim_dipole(axis,i1,j1-1,k1,i2,j2-1,k2);
    if (j2>1) 
      result += oo2zeta*(j2-1)*comp_prim_dipole(axis,i1,j1,k1,i2,j2-2,k2);
    if(axis==1) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2-1,k2);
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_dipole(axis,i1,j1,k1,i2,j2,k2-1);
    if (k1) 
      result += oo2zeta*i1*comp_prim_dipole(axis,i1,j1,k1-1,i2,j2,k2-1);
    if (k2>1) 
      result += oo2zeta*(k2-1)*comp_prim_dipole(axis,i1,j1,k1,i2,j2,k2-2);
    if(axis==2) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2,k2-1);
    return result;
    }

  return sMus[axis];
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
