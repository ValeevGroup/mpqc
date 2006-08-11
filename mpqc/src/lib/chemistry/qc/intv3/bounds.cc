//
// bounds.cc
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
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/int2e.h>

using namespace std;
using namespace sc;

#define COMPUTE_Q 1
#define COMPUTE_R 2

/* find the biggest number in the buffer */
static double
find_max(double *int_buffer,int nint)
{
  int i;
  double max = 0.0;
  for (i=0; i<nint; i++) {
    double val = int_buffer[i];
    if (val<0) val = -val;
    if (val > max) max = val;
    }
  return max;
  }

void
Int2eV3::int_init_bounds_nocomp()
{
  int i;
  int nshell=bs1_->nshell();
  int nsht=nshell*(nshell+1)/2;

  if (int_Qvec) free(int_Qvec);
  
  int_Qvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  used_storage_ += sizeof(int_bound_t)*nsht;
  if(int_Qvec==0) {
    ExEnv::errn() << scprintf("int_init_bounds_nocomp: cannot malloc int_Qvec: %d",
                     nsht)
         << endl;
    exit(1);
    }

  int_Rvec = 0;

  int_Q = int_bound_min;
  for (i=0; i<nsht; i++) int_Qvec[i] = 0;
}

void
Int2eV3::init_bounds()
{
  int_init_bounds_nocomp();
  compute_bounds(&int_Q,int_Qvec,COMPUTE_Q);
  }

void
Int2eV3::int_init_bounds_1der_nocomp()
{
  int i;
  int nshell=bs1_->nshell();
  int nsht=nshell*(nshell+1)/2;

  if (!int_derivative_bounds) {
    ExEnv::errn() << "requested der bounds but space not allocated" << endl;
    exit(1);
    }

  if (int_Qvec) free(int_Qvec);
  if (int_Rvec) free(int_Rvec);
  
  int_Qvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  int_Rvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  used_storage_ += sizeof(int_bound_t)*nsht*2;
  if((int_Qvec==0) || (int_Rvec==0)) {
    ExEnv::errn() << scprintf("int_init_bounds_1der_nocomp: cannot malloc int_{R,Q}vec: %d",nsht) << endl;
    exit(1);
    }

  int_Q = int_bound_min;
  int_R = int_bound_min;
  for (i=0; i<nsht; i++) int_Qvec[i] = int_Rvec[i] = 0;
  }

void
Int2eV3::int_bounds_comp(int s1, int s2)
{
  compute_bounds_shell(&int_Q,int_Qvec,COMPUTE_Q,s1,s2);
  }

void
Int2eV3::int_bounds_1der_comp(int s1, int s2)
{
  compute_bounds_shell(&int_Q,int_Qvec,COMPUTE_Q,s1,s2);
  compute_bounds_shell(&int_R,int_Rvec,COMPUTE_R,s1,s2);
  }

void
Int2eV3::init_bounds_1der()
{
  int_init_bounds_1der_nocomp();
  compute_bounds(&int_Q,int_Qvec,COMPUTE_Q);
  compute_bounds(&int_R,int_Rvec,COMPUTE_R);
  }

void
Int2eV3::done_bounds()
{
  if (int_Qvec) free(int_Qvec);
  int_Qvec = 0;
  }

void
Int2eV3::done_bounds_1der()
{
  if (int_Qvec) free(int_Qvec);
  if (int_Rvec) free(int_Rvec);
  int_Qvec = 0;
  int_Rvec = 0;
  }

int
Int2eV3::erep_4bound(int s1, int s2, int s3, int s4)
{
  if (!int_Qvec)
      return 256;
  
  int Qij;
  int Qkl;
  if (s1 >= 0 && s2 >= 0) {
      int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
      Qij = int_Qvec[ij];
    }
  else Qij = int_Q;
  if (s3 >=0 && s4 >= 0) {
      int kl=(s3>s4) ? ((s3*(s3+1))>>1)+s4 : ((s4*(s4+1))>>1)+s3;
      Qkl = int_Qvec[kl];
    }
  else Qkl = int_Q;

  return Qij+Qkl;
  }

int
Int2eV3::int_erep_2bound(int s1, int s2)
{
  if (!int_Qvec)
      return int_bound_max;
  
  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;

  return((int) int_Qvec[ij]);
  }

int
Int2eV3::int_erep_0bound_1der()
{
#if 0
  ExEnv::outn() << scprintf("int_erep_0bound_1der(): Q: %4d R: %4d\n", int_Q,int_R);
#endif
  return 1 + int_Q + int_R;
  }

int
Int2eV3::int_erep_2bound_1der(int s1, int s2)
{
  if (!int_Qvec || !int_Rvec)
      return int_bound_max;

  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
  int b1 = int_Qvec[ij] + int_R;
  int b2 = int_Q + int_Rvec[ij];

#if 0
  ExEnv::outn() << scprintf("int_erep_2bound_1der(%d,%d): Q: %4d R: %4d\n",s1,s2,
                   int_Qvec[ij],int_Rvec[ij]);
#endif

  /* The actual bound is Qij R + Q Rij
   * but since I'm using log base 2 I'll use
   * 2 * max (Qij R, Q Rij) -> 1 + max (Qij + R, Q + Rij)
   */

  return 1 + ((b1>b2)? b1 : b2);
  }

int
Int2eV3::erep_4bound_1der(int s1, int s2, int s3, int s4)
{
  if (!int_Qvec || !int_Rvec)
      return 256;

  int Qij, Qkl, Rij, Rkl;

  if (s1 >= 0 && s2 >= 0) {
      int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
      Qij = int_Qvec[ij];
      Rij = int_Rvec[ij];
    }
  else {
      Qij = int_Q;
      Rij = int_R;
    }
  if (s3 >= 0 && s4 >= 0) {
      int kl=(s3>s4) ? ((s3*(s3+1))>>1)+s4 : ((s4*(s4+1))>>1)+s3;
      Qkl = int_Qvec[kl];
      Rkl = int_Rvec[kl];
    }
  else {
      Qkl = int_Q;
      Rkl = int_R;
    }

  int b1 = Qij + Rkl;
  int b2 = Qkl + Rij;

#if 0
  ExEnv::outn() << scprintf("int_erep_4bound_1der(%d,%d,%d,%d): Q: %4d %4d R: %4d %4d\n",
                   s1,s2,s3,s4,
                   int_Qvec[ij],int_Qvec[kl],int_Rvec[ij],int_Rvec[kl]);
#endif

  /* The actual bound is Qij Rkl + Qkl Rij
   * but since I'm using log base 2 I'll use
   * 2 * max (Qij Rkl, Qkl Rij) -> 1 + max (Qij + Rkl, Qkl + Rij)
   */

  return 1 + ((b1>b2)? b1 : b2 );
  }

/* ripped off from clj's libintv2 */
/* (add subsequently ripped back on from ets's libdmtscf) */

/* Compute the partial bound arrays, either Q or R can be computed
 * with appropiate choice of flag. */
void
Int2eV3::compute_bounds(int_bound_t *overall, int_bound_t *vec, int flag)
{
  int sh1,sh2;

  if ((bs1_ != bs2_)||(bs1_ != bs3_)||(bs1_ != bs4_)) {
    ExEnv::errn() << scprintf("bounds.compute_bounds: all centers must be the same")
         << endl;
    exit(1);
    }

  int nshell=bs1_->nshell();
  int nsht=(nshell*(nshell+1))/2;

  int me = grp_->me();
  int n = grp_->n();

  for (int i=0; i<nsht; i++) vec[i] = 0;

  *overall = int_bound_min;
  int sh12 = 0;
  for(sh1=0; sh1 < bs1_->nshell() ; sh1++) {
    for(sh2=0; sh2 <= sh1 ; sh2++,sh12++) {
      if (sh12%n == me) compute_bounds_shell(overall,vec,flag,sh1,sh2);
      }
    }

  grp_->sum(vec,nsht);
  grp_->max(overall,1);
  }

/* Compute the partial bound arrays, either Q or R can be computed
 * with appropiate choice of flag. */
void
Int2eV3::compute_bounds_shell(int_bound_t *overall, int_bound_t *vec,
                              int flag, int sh1, int sh2)
{
  int nint;
  int shellij;
  int shells[4],size[4];
  double max;
  double tol = pow(2.0,double(int_bound_min));
  double loginv = 1.0/log(2.0);
  int old_int_integral_storage = int_integral_storage;
  int_integral_storage = 0;

  int old_perm = permute();
  set_permute(0);
  int old_red = redundant();
  set_redundant(1);

  if ((bs1_ != bs2_)||(bs1_ != bs3_)||(bs1_ != bs4_)) {
    ExEnv::errn() << scprintf("bounds.compute_bounds: all centers must be the same")
         << endl;
    exit(1);
    }

  if (sh1<sh2) {
    int tmp = sh1;
    sh1 = sh2;
    sh2 = tmp;
    }

  shellij= ((sh1*(sh1+1))>>1) + sh2;
    shells[0]=shells[2]=sh1;
      shells[1]=shells[3]=sh2;

      if (flag == COMPUTE_Q) {
        erep(shells,size);
        nint = size[0]*size[1]*size[0]*size[1];
        max = find_max(int_buffer,nint);
#if 0
        ExEnv::outn() << scprintf("max for %d %d (size %d) is %15.11f\n", sh1, sh2, nint, max);
#endif
        }
      else if (flag == COMPUTE_R) {
        double max1,max2;
        int_erep_bound1der(0,sh1,sh2,&nint);
        max1 = find_max(int_buffer,nint);
#if 0
        ExEnv::outn() << scprintf("bound(%d) for (%d,%d) is %12.8f int_buffer =",
                         flag,sh1,sh2,max1);
        for (i=0; (i<nint)&&(i<27); i++)
          ExEnv::outn() << scprintf(" %12.8f",int_buffer[i]);
        if (nint > 27) ExEnv::outn() << scprintf(" ...");
        ExEnv::outn() << scprintf("\n");
#endif
        int_erep_bound1der(0,sh2,sh1,&nint);
        max2 = find_max(int_buffer,nint);
        max = (max1>max2)?max1:max2;
        }
      else {
        ExEnv::outn() << scprintf("bad bound flag\n"); exit(1);
        }

    /* Compute the partial bound value. */
      max = sqrt(max);
      if (max>tol) {
        vec[shellij] = (int_bound_t) ceil(log(max)*loginv);
        }
      else {
        vec[shellij] = (int_bound_t) int_bound_min;
        }

    /* Multiply R contributions by a factor of two to account for
     * fact that contributions from both centers must be accounted
     * for. */
      if (flag == COMPUTE_R) vec[shellij]++;
      if (vec[shellij]>*overall) *overall = vec[shellij];
#if 0
      ExEnv::outn() << scprintf("bound(%d) for (%d,%d) is %4d int_buffer =",
                       flag,sh1,sh2,vec[shellij]);
      for (i=0; (i<nint)&&(i<27); i++) ExEnv::outn() << scprintf(" %12.8f",int_buffer[i]);
      if (nint > 27) ExEnv::outn() << scprintf(" ...");
      ExEnv::outn() << scprintf("\n");
#endif
  int_integral_storage = old_int_integral_storage;

  set_permute(old_perm);
  set_redundant(old_red);

  }

/* This function is used to convert a double to its log base 2 rep
 * for use in bound computations. */
int
Int2eV3::bound_to_logbound(double value)
{
  double tol = pow(2.0,double(int_bound_min));
  double loginv = 1.0/log(2.0);
  int_bound_t res;

  if (value > tol) res = (int_bound_t) ceil(log(value)*loginv);
  else res = int_bound_min;
  return res;
  }

/* This function is used to convert a double from its log base 2 rep. */
double
Int2eV3::logbound_to_bound(int value)
{
  return pow(2.0,(double)value);
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
