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
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/fjt.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/tformv3.h>

#define IN(i,j) ((i)==(j)?1:0)
#define SELECT(x1,x2,x3,s) (((s)==0)?x1:(((s)==1)?(x2):(x3)))

/* ------------ Initialization of 1e routines. ------------------- */
/* This routine returns a buffer large enough to hold a shell doublet
 * of integrals (if order == 0) or derivative integrals (if order == 1).
 */
void
Int1eV3::int_initialize_1e(int flags, int order)
{
  int jmax1,jmax2,jmax;
  int scratchsize,nshell2;

  /* The efield routines look like derivatives so bump up order if
   * it is zero to allow efield integrals to be computed.
   */
  if (order == 0) order = 1;

  jmax1 = bs1_->max_angular_momentum();
  jmax2 = bs2_->max_angular_momentum();
  jmax = jmax1 + jmax2;

  fjt_ = new FJT(jmax + 2*order);

  nshell2 = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell();

  if (order == 0) {
    init_order = 0;
    scratchsize = nshell2;
    }
  else if (order == 1) {
    init_order = 1;
    scratchsize = nshell2*3;
    }
  else {
    cerr << scprintf("int_initialize_1e: invalid order: %d\n",order);
    exit(1);
    }

  buff = (double *) malloc(scratchsize*sizeof(double));
  cartesianbuffer = (double *) malloc(scratchsize*sizeof(double));

  }

void
Int1eV3::int_done_1e()
{
  init_order = -1;
  free(buff);
  free(cartesianbuffer);
  buff = 0;
  cartesianbuffer = 0;
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
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

  intv3_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the overlap ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::overlap_1der(int ish, int jsh,
                      int idercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_overlap: one electron routines are not init'ed\n");
    exit(1);
    }

  RefGaussianBasisSet dercs;
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
  cout << scprintf("zeroing %d*%d*3 elements of buff\n",ni,nj);
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
                                      RefGaussianBasisSet dercs, int centernum)
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
      ss =   pow(3.141592653589793/(exp1+exp2),1.5)
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
  ss =   pow(3.141592653589793/(gshell1->exponent(prim1)
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
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

  intv3_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
  }

void
Int1eV3::int_accum_shell_kinetic(int ish, int jsh)
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
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
  intv3_accum_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the kinetic energy derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
void
Int1eV3::int_accum_shell_kinetic_1der(int ish, int jsh,
                                      RefGaussianBasisSet dercs, int centernum)
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
                          RefGaussianBasisSet dercs, int centernum,
                          double (Int1eV3::*shell_function)
                          (int,int,int,int,int,int,int,int))
{
  int i;
  int gc1,gc2;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
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

  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, gshell1, gshell2);
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
      ss =   pow(3.141592653589793/(gshell1->exponent(i)
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
                                      RefGaussianBasisSet dercs, int centernum)
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
                                          RefGaussianBasisSet dercs,
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
                                         RefGaussianBasisSet dercs,
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
    accum_shell_efield(buff,ish,jsh);
    scale_shell_result = 0;
    }
  else if (bs2_ == dercs) {
    scale_shell_result = 1;
    result_scale_factor= -bs2_->molecule()->charge(centernum);
    for (int xyz=0; xyz<3; xyz++) {
        C[xyz] = bs2_->r(centernum,xyz);
      }
    accum_shell_efield(buff,ish,jsh);
    scale_shell_result = 0;
    }

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the non Hellman-Feynman part is computed.
 */
void
Int1eV3::int_accum_shell_nuclear_nonhf_1der(int ish, int jsh,
                                            RefGaussianBasisSet dercs,
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
    accum_shell_1der(buff,ish,jsh,dercs,centernum,
                     &Int1eV3::comp_shell_nuclear);
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
      accum_shell_1der(buff,ish,jsh,dercs,centernum,
                       &Int1eV3::comp_shell_nuclear);
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    cerr << scprintf("accum_shell_efield: one electron routines are not init'ed\n");
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

  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the efield integrals between functions in two shells.
 * The result is placed in the buffer in the form bf1 x y z, bf2
 * x y z, etc.
 */
void
Int1eV3::efield(int ish, int jsh, double *position)
{
  scale_shell_result = 0;
  int xyz;

  for (xyz=0; xyz<3; xyz++) {
      C[xyz] = position[xyz];
    }

  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    cerr << scprintf("Int1eV3::efield one electron routines are not ready\n");
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

  intv3_transform_1e_xyz(cartesianbuffer, buff, gshell1, gshell2);
}

/* This computes the nuc rep energy integrals between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::nuclear(int ish, int jsh)
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_nuclear: one electron routines are not init'ed\n");
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

  intv3_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double tmp;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_pointcharge: one electron routines are not init'ed\n");
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

  intv3_accum_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;

  if (!(init_order >= 0)) {
    cerr << scprintf("Int1eV3::point_charge: one electron routines are not init'ed\n");
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

  intv3_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
  }


/* This computes the 1e Hamiltonian integrals between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::hcore(int ish, int jsh)
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int cart1,cart2;
  int gc1,gc2;

  if (!(init_order >= 0)) {
    cerr << scprintf("hcore: one electron routines are not init'ed\n");
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

  intv3_transform_1e(cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the 1e Hamiltonian deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
void
Int1eV3::hcore_1der(int ish, int jsh,
                    int idercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  RefGaussianBasisSet dercs;
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
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_kinetic: one electron routines are not init'ed\n");
    exit(1);
    }

  RefGaussianBasisSet dercs;
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
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  RefGaussianBasisSet dercs;
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
                                   RefGaussianBasisSet dercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_nuclear_hf_1der: one electron routines are not init'ed\n");
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
                                      RefGaussianBasisSet dercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    cerr << scprintf("int_shell_nuclear_nonhf_1der: one electron routines are not init'ed\n");
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
  cout << scprintf("int_shell_nuclear_nonhf_1der: zeroing %d doubles in buff\n",ni*nj*3);
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
      auxcoef =   2.0 * 3.141592653589793/(gshell1->exponent(i)
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
      auxcoef =   2.0 * 3.141592653589793/(gshell1->exponent(i)
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
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int im,jm,km;
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
  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, gshell1, gshell2);
  }

/* This computes the dipole integrals between functions in two shells.
 * The result is placed in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The last arg, com, is the origin of the coordinate
 * system used to compute the dipole moment.
 */
void
Int1eV3::dipole(int ish, int jsh, double *com)
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int im,jm,km;
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
  intv3_transform_1e_xyz(cartesianbuffer, buff, gshell1, gshell2);
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
      ss =   pow(3.141592653589793/(exp1+exp2),1.5)
           * exp(- oozeta * exp1 * exp2 * AmB2);
      sMus = ss * PmC[mu];
      tmp     =  gshell1->coefficient_unnorm(gc1,i)
               * gshell2->coefficient_unnorm(gc2,j);
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      dipole[0] += tmp * comp_prim_dipole(1,0,0,i1,j1,k1,i2,j2,k2);
      dipole[1] += tmp * comp_prim_dipole(0,1,0,i1,j1,k1,i2,j2,k2);
      dipole[2] += tmp * comp_prim_dipole(0,0,1,i1,j1,k1,i2,j2,k2);
      }
    }

  }

double
Int1eV3::comp_prim_dipole(int im, int jm, int km,
                          int i1, int j1, int k1,
                          int i2, int j2, int k2)
{
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2,j2,k2);
    if (i2) 
      result += oo2zeta*i2*comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2-1,j2,k2);
    if (i1>1)
      result += oo2zeta*(i1-1)*comp_prim_dipole(im,jm,km,i1-2,j1,k1,i2,j2,k2);
    if(im) result += oo2zeta*comp_prim_overlap(i1-1,j1,k1,i2,j2,k2);
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2,k2);
    if (j2) 
      result += oo2zeta*j2*comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2-1,k2);
    if (j1>1) 
      result += oo2zeta*(j1-1)*comp_prim_dipole(im,jm,km,i1,j1-2,k1,i2,j2,k2);
    if(jm) result += oo2zeta*comp_prim_overlap(i1,j1-1,k1,i2,j2,k2);
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2);
    if (k2) 
      result += oo2zeta*k2*comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2-1);
    if (k1>1) 
      result += oo2zeta*(k1-1)*comp_prim_dipole(im,jm,km,i1,j1,k1-2,i2,j2,k2);
    if(km) result += oo2zeta*comp_prim_overlap(i1,j1,k1-1,i2,j2,k2);
    return result;
    }
  if (i2) {
    result = PmB[0] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2-1,j2,k2);
    if (i1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2-1,j2,k2);
    if (i2>1) 
      result += oo2zeta*(i2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2-2,j2,k2);
    if(im) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2-1,k2);
    if (j1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2-1,k2);
    if (j2>1) 
      result += oo2zeta*(j2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2-2,k2);
    if(jm) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2-1,k2);
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2,k2-1);
    if (k1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2-1);
    if (k2>1) 
      result += oo2zeta*(k2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2,k2-2);
    if(km) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2,k2-1);
    return result;
    }

  return sMus;
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
