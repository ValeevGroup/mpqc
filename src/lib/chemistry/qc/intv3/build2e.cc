//
// build2e.cc
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
#include <cassert>
#include <math.h>

#include <mpqc_config.h>
#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/int2e.h>

using namespace std;
using namespace sc;

#define CHECK_STACK_ALIGNMENT 0
#if CHECK_STACK_ALIGNMENT
static void
stack_alignment_error(void *ptr, const char *where)
{
  ExEnv::outn() << "UNALIGNED STACK: " << where << ": " << ptr << endl;
}
static inline void
stack_alignment_check(void *ptr, const char *where)
{
  if ((unsigned)ptr & 7) stack_alignment_error(ptr,where);
}
#else
#  define stack_alignment_check(ptr,where)
#endif

  /* MG is the maximum angular momentum for which we will use
   * the generated build routines. It is defined in oint3/build.h */
#define MINA(x) (((x)<MG)?(x):MG)

static inline void 
iswtch(int *i, int *j)
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

/* This initializes the build routines.  It is called from
 * int_initialize_erep.  This allocates storage for the
 * intermediate integrals. */
void
Int2eV3::int_init_buildgc(int order,
                          int am1, int am2, int am3, int am4,
                          int nc1, int nc2, int nc3, int nc4)
{
  int *jmax_for_con;
  int am12;
  int am34;
  int am;
  int i,j,k,l,m;
  int ci,cj,ck,cl;
  int e0f0_con_int_bufsize;
  double *e0f0_con_int_buf;
  int int_v_bufsize, int_v0_bufsize;
  double *int_v_buf, *int_v0_buf;

  used_storage_build_ = 0;

  /* Convert the am1-4 to their canonical ordering. */
  if (am2>am1) {
    iswtch(&am1,&am2);
    iswtch(&nc1,&nc2);
    }
  if (am4>am3) {
    iswtch(&am3,&am4);
    iswtch(&nc3,&nc4);
    }
  if ((am3 > am1)||((am3 == am1)&&(am4 > am2))) {
    iswtch(&am1,&am3);
    iswtch(&nc1,&nc3);
    iswtch(&am2,&am4);
    iswtch(&nc2,&nc4);
    }

  /* If the center permutation 1<->3 and 2<->4 is performed, then
   * we may need the am for center 2 to be as big as for center 4. */
  if (am4 > am2) am2 = am4;

  /* As far as this routine knows the biggest nc can end up anywhere. */
  if (nc2>nc1) nc1 = nc2;
  if (nc3>nc1) nc1 = nc3;
  if (nc4>nc1) nc1 = nc4;
  nc2 = nc3 = nc4 = nc1;

  jmax_for_con = (int *) malloc(sizeof(int)*nc1);
  // storage for jmax_for_con is not counted since it is freed below
  for (i=0; i<nc1; i++) {
    int tmp;
    jmax_for_con[i] = bs1_->max_am_for_contraction(i);
    if (  (bs2_ != bs1_)
        &&((tmp=(int_unit2?0:bs2_->max_am_for_contraction(i)))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (bs3_ != bs1_) && (bs3_ != bs2_)
        &&((tmp=bs3_->max_am_for_contraction(i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (bs4_ != bs1_) && (bs4_ != bs2_) && (bs4_ != bs3_)
        &&((tmp=(int_unit4?0:bs4_->max_am_for_contraction(i)))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    }

  /* If derivatives are needed, then am1 can be bigger. */
  if (order==1) am1++;
  /* To compute derivative integral bounds, am3 can be bigger also. */
  if (order==1 && int_derivative_bounds) am3++;

  am12 = am1 + am2;
  am34 = am3 + am4;
  am = am12 + am34;

  /* Allocate the intlist. */
  contract_length.set_dim(am12+1,am34+1,am34+1);
  build.int_v_list.set_dim(am12+1,am34+1,am+1);
  used_storage_build_ += contract_length.nbyte();
  used_storage_build_ += build.int_v_list.nbyte();
#if CHECK_INTEGRAL_ALGORITHM
  ExEnv::outn() << "contract_length: " << contract_length.nbyte() << endl;
  ExEnv::outn() << "int_v_list: " << build.int_v_list.nbyte() << endl;
#endif

  /* Set all slots to 0 */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34; k++) {
        build.int_v_list(i,j,k) = 0;
        }
      }
    }

  for (i=0; i<=am12; i++) {
      for (j=0; j<=am34; j++) {
          for (k=0; k<=am34; k++) {
              contract_length(i,j,k) = 0;
              for (l=j; l<=k; l++) {
                  contract_length(i,j,k) += INT_NCART(i)*INT_NCART(l);
                }
            }
        }
    }

  /* Compute the size of the buffer for the primitive integrals. */
  int_v_bufsize = 0;
  int_v0_bufsize = 0;
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      int_v0_bufsize += INT_NCART(i)*INT_NCART(j);
      for (k=0; k<=am12+am34-i-j; k++) {
        int_v_bufsize += INT_NCART(i)*INT_NCART(j);
        }
      }
    }

  int_v0_buf = (double*) allocate<char>(sizeof(double)*int_v_bufsize);
  used_storage_build_ += sizeof(double)*int_v_bufsize;
  if (!int_v0_buf) {
    ExEnv::errn() << scprintf("couldn't allocate all integral intermediates\n");
    fail();
    }
  add_store(int_v0_buf);
  int_v_buf = &int_v0_buf[int_v0_bufsize];

  /* Allocate storage for the needed slots. */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      build.int_v_list(i,j,0) = int_v0_buf;
      int_v0_buf += INT_NCART(i)*INT_NCART(j);
      for (k=1; k<=am12+am34-i-j; k++) {
        build.int_v_list(i,j,k) = int_v_buf;
        int_v_buf += INT_NCART(i)*INT_NCART(j);
        }
      }
    }


  /* Allocate storage for the contracted integrals (these are the output
   * of the build routines). */
  /* The ci, etc, indices refer to which set of contraction
   * coefficients we are using. */
  e0f0_con_int_bufsize = 0;
  e0f0_con_ints_array = new IntV3Arraydoublep2***[nc1];
  used_storage_build_ += sizeof(IntV3Arraydoublep2***)*nc1;
  for (ci=0; ci<nc1; ci++) {
    e0f0_con_ints_array[ci] = new IntV3Arraydoublep2**[nc2];
    used_storage_build_ += sizeof(IntV3Arraydoublep2**)*nc2;
    for (cj=0; cj<nc2; cj++) {
      e0f0_con_ints_array[ci][cj] = new IntV3Arraydoublep2*[nc3];
      used_storage_build_ += sizeof(IntV3Arraydoublep2*)*nc3;
      for (ck=0; ck<nc3; ck++) {
        e0f0_con_ints_array[ci][cj][ck] = new IntV3Arraydoublep2[nc4];
        used_storage_build_ += sizeof(IntV3Arraydoublep2)*nc4;
        for (cl=0; cl<nc4; cl++) {
  int am12_for_con;
  int am34_for_con;

  am12_for_con = jmax_for_con[ci] + jmax_for_con[cj] + order;
  if ((jmax_for_con[ck]!=am3)||(jmax_for_con[cl]!=am4)) {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl] + order;
    }
  else {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl];
    }

#if CHECK_INTEGRAL_ALGORITHM
  ExEnv::outn() << "am12_for_con: " << am12_for_con << endl;
  ExEnv::outn() << "am34_for_con: " << am34_for_con << endl;
#endif

  e0f0_con_ints_array[ci][cj][ck][cl].set_dim(am12_for_con+1,am34_for_con+1);
  used_storage_build_ += e0f0_con_ints_array[ci][cj][ck][cl].nbyte();
#if CHECK_INTEGRAL_ALGORITHM
  ExEnv::outn() << "e0f0_con_ints_array: "
       << e0f0_con_ints_array[ci][cj][ck][cl].nbyte()
       << endl;
#endif

  /* Count how much storage for the integrals is needed. */
  for (i=0; i<=am12_for_con; i++) {
    for (k=0; k<=am34_for_con; k++) {
      int s =  INT_NCART(i)
               * INT_NCART(k);
      e0f0_con_int_bufsize += s;
      }
    }
          }
        }
      }
    }
  e0f0_con_int_buf = (double*)allocate<char>(sizeof(double)*e0f0_con_int_bufsize);
  used_storage_build_ += e0f0_con_int_bufsize * sizeof(double);
#if CHECK_INTEGRAL_ALGORITHM
  ExEnv::outn() << "e0f0_int_buf: " << e0f0_con_int_bufsize * sizeof(double) << endl;
#endif
  if (!e0f0_con_int_buf) {
    ExEnv::errn() << scprintf("couldn't allocate contracted integral storage\n");
    fail();
    }
  add_store(e0f0_con_int_buf);
  /* Allocate storage for the integrals which will be used by the shift
   * routine. */
  for (ci=0; ci<nc1; ci++) {
    for (cj=0; cj<nc2; cj++) {
      for (ck=0; ck<nc3; ck++) {
        for (cl=0; cl<nc4; cl++) {
  int am12_for_con;
  int am34_for_con;

  am12_for_con = jmax_for_con[ci] + jmax_for_con[cj] + order;
  if ((jmax_for_con[ck]!=am3)||(jmax_for_con[cl]!=am4)) {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl] + order;
    }
  else {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl];
    }

  for (i=0; i<=am12_for_con; i++) {
    for (k=0; k<=am34_for_con; k++) {
      e0f0_con_ints_array[ci][cj][ck][cl](i,k) = 0;
      }
    }
  for (i=0; i<=am12_for_con; i++) {
    for (k=0; k<=am34_for_con; k++) {
/* If there are Pople style s=p shells and the shells are ordered
 * first s and then p and there are no p or d shells on the molecule,
 * then this algorithm would will allocate a little more storage
 * than needed.  General contraction should be ordered high to
 * low angular momentum for this reason. */
      e0f0_con_ints_array[ci][cj][ck][cl](i,k)
        = e0f0_con_int_buf;
      e0f0_con_int_buf +=  INT_NCART(i)
                           * INT_NCART(k);
      }
    }
          }
        }
      }
    }

  /* Initialize the build_routine array. */
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      for (k=0; k<4; k++) {
        for (l=0; l<4; l++) {
          for (m=0; m<2; m++) {
    build_routine[i][j][k][l][m] = &BuildIntV3::impossible_integral;
            }
          }
        }
      }
    }

#define ASSIGN_BUILD(ii,j,k,l) \
  build_routine[ii][j][k][l][0]= &BuildIntV3::i ## ii ## j ## k ## l ;\
  build_routine[ii][j][k][l][1]= &BuildIntV3::i ## ii ## j ## k ## l ## eAB;
#if (MG == 1) || (MG == 2) || (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,1,0,0)
  ASSIGN_BUILD(0,1,0,1)
  ASSIGN_BUILD(0,1,1,1)
  ASSIGN_BUILD(1,1,0,0)
  ASSIGN_BUILD(1,1,1,1)
#endif

#if (MG == 2) || (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,2,0,0)
  ASSIGN_BUILD(0,2,0,1)
  ASSIGN_BUILD(0,2,0,2)
  ASSIGN_BUILD(0,2,1,1)
  ASSIGN_BUILD(0,2,1,2)
  ASSIGN_BUILD(0,2,2,2)
  ASSIGN_BUILD(1,2,0,0)
  ASSIGN_BUILD(1,2,0,1)
  ASSIGN_BUILD(1,2,1,1)
  ASSIGN_BUILD(1,2,1,2)
  ASSIGN_BUILD(1,2,2,2)
  ASSIGN_BUILD(2,2,0,0)
  ASSIGN_BUILD(2,2,0,1)
  ASSIGN_BUILD(2,2,1,1)
  ASSIGN_BUILD(2,2,2,2)
#endif

#if (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,3,0,0)
  ASSIGN_BUILD(0,3,0,1)
  ASSIGN_BUILD(0,3,0,2)
  ASSIGN_BUILD(0,3,0,3)
  ASSIGN_BUILD(0,3,1,1)
  ASSIGN_BUILD(0,3,1,2)
  ASSIGN_BUILD(0,3,1,3)
  ASSIGN_BUILD(0,3,2,2)
  ASSIGN_BUILD(0,3,2,3)
  ASSIGN_BUILD(0,3,3,3)
  ASSIGN_BUILD(1,3,0,0)
  ASSIGN_BUILD(1,3,0,1)
  ASSIGN_BUILD(1,3,0,2)
  ASSIGN_BUILD(1,3,1,1)
  ASSIGN_BUILD(1,3,1,2)
  ASSIGN_BUILD(1,3,1,3)
  ASSIGN_BUILD(1,3,2,2)
  ASSIGN_BUILD(1,3,2,3)
  ASSIGN_BUILD(1,3,3,3)
  ASSIGN_BUILD(2,3,0,0)
  ASSIGN_BUILD(2,3,0,1)
  ASSIGN_BUILD(2,3,0,2)
  ASSIGN_BUILD(2,3,1,1)
  ASSIGN_BUILD(2,3,1,2)
  ASSIGN_BUILD(2,3,2,2)
  ASSIGN_BUILD(2,3,2,3)
  ASSIGN_BUILD(2,3,3,3)
  ASSIGN_BUILD(3,3,0,0)
  ASSIGN_BUILD(3,3,0,1)
  ASSIGN_BUILD(3,3,0,2)
  ASSIGN_BUILD(3,3,1,1)
  ASSIGN_BUILD(3,3,1,2)
  ASSIGN_BUILD(3,3,2,2)
  ASSIGN_BUILD(3,3,3,3)
#endif

  free(jmax_for_con);
  saved_am12 = am12;
  saved_am34 = am34;
  saved_ncon = nc1;

  used_storage_ += used_storage_build_;
#if CHECK_INTEGRAL_ALGORITHM
  ExEnv::outn() << "used_storage_build: " << used_storage_build_ << endl;
#endif
  }

void
Int2eV3::int_done_buildgc()
{
  int ci,cj,ck;

  used_storage_ -= used_storage_build_;
  used_storage_build_ = 0;

  free_store();

  for (ci=0; ci<saved_ncon; ci++) {
    for (cj=0; cj<saved_ncon; cj++) {
      for (ck=0; ck<saved_ncon; ck++) {
        delete[] e0f0_con_ints_array[ci][cj][ck];
        }
      delete[] e0f0_con_ints_array[ci][cj];
      }
    delete[] e0f0_con_ints_array[ci];
    }
  delete[] e0f0_con_ints_array;

  }

/* add_store maintains a list of free storage allocated by int_init_buildgc */
void
Int2eV3::add_store(void *p)
{
  if (!store) {
    store = (store_list_t*) malloc(sizeof(store_list_t));
    MPQC_ASSERT(store);
    store->p = 0;
    n_store_last = 0;
    }
  if (n_store_last >= STORAGE_CHUNK) {
    store_list_t* tmp = (store_list_t*) malloc(sizeof(store_list_t));
    MPQC_ASSERT(tmp);
    tmp->p = store;
    store = tmp;
    n_store_last = 0;
    }
  store->data[n_store_last++] = p;
  }

/* free_store frees the memory that add_store keeps track of */
void
Int2eV3::free_store()
{
  _free_store(store,n_store_last);
  store = 0;
  }

void
Int2eV3::_free_store(store_list_t* s, int n)
{
  int i;
  if (!s) return;
  for (i=0; i<n; i++) {
    deallocate(static_cast<char*>(s->data[i]));
    }
  _free_store(s->p,STORAGE_CHUNK);
  free(s);
  }


void
Int2eV3::int_buildgcam(int minam1, int minam2, int minam3, int minam4,
                       int maxam1, int maxam2, int maxam3, int maxam4,
                       int dam1, int dam2, int dam3, int dam4,
                       int sh1, int sh2, int sh3, int sh4,
                       int eAB)
{
  int k,m,n;
  int ci,cj,ck,cl;
  int maxam12,maxam34;
  int nc1,nc2,nc3,nc4;

  if (maxam1<0 || maxam2<0 || maxam3<0 || maxam4<0) return;
  if (minam1<0) minam1=0;
  if (minam2<0) minam2=0;
  if (minam3<0) minam3=0;
  if (minam4<0) minam4=0;

  maxam12 = maxam1 + maxam2;
  maxam34 = maxam3 + maxam4;

  nc1 = pbs1_->shell(sh1).ncontraction();
  if (pbs2_ == 0) nc2 = 1;
  else nc2 = pbs2_->shell(sh2).ncontraction();
  nc3 = pbs3_->shell(sh3).ncontraction();
  if (pbs4_ == 0) nc4 = 1;
  else nc4 = pbs4_->shell(sh4).ncontraction();

  /* Zero the target contracted integrals that the build routine
   * will accumulate into. */
  for (m=minam1; m<=maxam12; m++) {
    for (n=minam3; n<=maxam34; n++) {
  int nm_cart = INT_NCART(m)*INT_NCART(n);
  for (ci=0; ci<nc1; ci++) {
    if (m < int_shell1->am(ci)+dam1) continue;
    for (cj=0; cj<nc2; cj++) {
      if (int_shell1->am(ci)+dam1+int_shell2->am(cj)+dam2 < m)
        continue;
      for (ck=0; ck<nc3; ck++) {
        if (n < int_shell3->am(ck)+dam3) continue;
        for (cl=0; cl<nc4; cl++) {
          if (int_shell3->am(ck)+dam3 +int_shell4->am(cl)+dam4 < n)
            continue;
          double *tmp = e0f0_con_ints_array[ci][cj][ck][cl](m,n);
          for (int ii=0; ii<nm_cart; ii++) tmp[ii] = 0.0;
          }
        }
      }
    }
      }
    }

  gen_shell_intermediates(sh1,sh2,sh3,sh4);

  /* If enough of the functions come from generalized contractions
   * to make it workwhile, then don't do redundant primitives
   * at the additional cost of slower normalization computations.
   */
  if (nc1 + nc2 + nc3 + nc4 > 4)
    build_using_gcs(nc1,nc2,nc3,nc4,
                    minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB);
  else
    build_not_using_gcs(nc1,nc2,nc3,nc4,
                    minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB);
  }

void
Int2eV3::build_not_using_gcs(int nc1, int nc2, int nc3, int nc4,
                             int minam1, int minam3, int maxam12, int maxam34,
                             int dam1, int dam2, int dam3, int dam4, int eAB)
{
  int i,j,k,l,m;
  int ci,cj,ck,cl;
  double *bufferprim;

#if 0
  ExEnv::outn() << scprintf("not_gcs: %d%d%d%d\n",
                   int_expweight1,
                   int_expweight2,
                   int_expweight3,
                   int_expweight4
      );
#endif

          /* Sum thru all possible contractions. */
  for (ci=0; ci<nc1; ci++) {
    int mlower = int_shell1->am(ci) + dam1;
    if (mlower < 0) continue;
    IntV3Arraydoublep2 ***e0f0_i = e0f0_con_ints_array[ci];
    for (cj=0; cj<nc2; cj++) {
      int mupper = mlower + int_shell2->am(cj) + dam2;
      if (mupper < mlower) continue;
      if (mlower < minam1) mlower = minam1;
      if (mupper > maxam12) mupper = maxam12;
      IntV3Arraydoublep2 **e0f0_ij = e0f0_i[cj];
      for (ck=0; ck<nc3; ck++) {
        int nlower = int_shell3->am(ck) + dam3;
        if (nlower < 0) continue;
        IntV3Arraydoublep2 *e0f0_ijk = e0f0_ij[ck];
        for (cl=0; cl<nc4; cl++) {
          int nupper = nlower + int_shell4->am(cl) + dam4;
          if (nupper < nlower) continue;
          if (nlower < minam3) nlower = minam3;
          if (nupper > maxam34) nupper = maxam34;

  /* Loop over the primitives. */
  for (i=0; i<int_shell1->nprimitive(); i++) {
    double coef0;
    coef0 = int_shell1->coefficient_unnorm(ci,i);
    if (int_expweight1) coef0 = coef0
                                    * int_shell1->exponent(i);
    /* This factor of two comes from the derivative integral formula. */
    if (int_expweight1) coef0 *= 2.0;
    if (int_expweight2) coef0 *= 2.0;
    if (int_expweight3) coef0 *= 2.0;
    if (int_expweight4) coef0 *= 2.0;
    if (int_store1) opr1 = int_shell_to_prim[osh1] + i;
    for (j=0; j<int_shell2->nprimitive(); j++) {
      double coef1;
      coef1 = int_shell2->coefficient_unnorm(cj,j);
      if (int_expweight2) coef1 *=  coef0
                                      * int_shell2->exponent(j);
      else                     coef1 *= coef0;
      if (int_store1) opr2 = int_shell_to_prim[osh2] + j;
      for (k=0; k<int_shell3->nprimitive(); k++) {
        double coef2;
        coef2 = int_shell3->coefficient_unnorm(ck,k);
        if (int_expweight3) coef2 *=  coef1
                                        * int_shell3->exponent(k);
        else                     coef2 *= coef1;
        if (int_store1) opr3 = int_shell_to_prim[osh3] + k;
        for (l=0; l<int_shell4->nprimitive(); l++) {
          double coef3;
          coef3 = int_shell4->coefficient_unnorm(cl,l);
          if (int_expweight4) coef3 *=  coef2
                                          * int_shell4->exponent(l);
          else                     coef3 *= coef2;
          if (int_store1) opr4 = int_shell_to_prim[osh4] + l;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates_with_norm(i,j,k,l, maxam12+maxam34,coef3);

          /* Generate the target integrals. */
          if ((maxam12 == 0) && (maxam34 == 0)) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            if (build_routine[minam1]
                             [maxam12]
                             [minam3]
                             [maxam34][eAB]==&BuildIntV3::impossible_integral){
              ExEnv::errn() << scprintf("trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if (!(build.*build_routine[minam1]
                                      [maxam12]
                                      [minam3]
                                      [maxam34][eAB])()) {
              ExEnv::outn() << "build2e.cc: did not succeed in building all integrals"
                   << endl;
              abort();
              }
            }
          else {
            blockbuildprim(minam1,maxam12,minam3,maxam34);
            }

          /* Contract the primitive target integrals. */
          /* Throw out all unneeded contractions. */
          if (i||j||k||l) {
            for (m=mlower; m<=mupper; m++) {
              int o;
              int sizec = contract_length(m,nlower,nupper);
              double *RESTRICT con_ints = e0f0_ijk[cl](m,nlower);
              bufferprim = build.int_v_list(m,nlower,0);

              for (o=sizec; o!=0; o--) {
                *con_ints++ += *bufferprim++;
                }

              }
            }
          else {
            // for the first primitive write to con_ints rather
            // than accumulate into it
            for (m=mlower; m<=mupper; m++) {
              int o;
              int sizec = contract_length(m,nlower,nupper);
              double *RESTRICT con_ints = e0f0_ijk[cl](m,nlower);
              bufferprim = build.int_v_list(m,nlower,0);

              for (o=sizec; o!=0; o--) {
                *con_ints++ = *bufferprim++;
                }

              }
            }

          }
        }
      }
    }

          }
        }
      }
    }

  }

void
Int2eV3::build_using_gcs(int nc1, int nc2, int nc3, int nc4,
                         int minam1, int minam3, int maxam12, int maxam34,
                         int dam1, int dam2, int dam3, int dam4, int eAB)
{
  int i,j,k,l,m;
  int ci,cj,ck,cl;
  int maxam1234=maxam12+maxam34;
  double coef0,coef1,coef2,coef3;
  double ishl1expi=1.0, ishl2expj=1.0, ishl3expk=1.0;
  double *bufferprim;
  double c0scale;

  /* Loop over the primitives. */
  for (i=0; i<int_shell1->nprimitive(); i++) {
    if (int_store1) opr1 = int_shell_to_prim[osh1] + i;
    if (int_expweight1) ishl1expi=2.0*int_shell1->exponent(i);

    for (j=0; j<int_shell2->nprimitive(); j++) {
      if (int_store1) opr2 = int_shell_to_prim[osh2] + j;
      ishl2expj = (int_expweight2) ? 
                        2.0*int_shell2->exponent(j)*ishl1expi : ishl1expi;

      for (k=0; k<int_shell3->nprimitive(); k++) {
        if (int_store1) opr3 = int_shell_to_prim[osh3] + k;
        ishl3expk = (int_expweight3) ? 
                        2.0*int_shell3->exponent(k)*ishl2expj : ishl2expj;

        for (l=0; l<int_shell4->nprimitive(); l++) {
          if (int_store1) opr4 = int_shell_to_prim[osh4] + l;
          c0scale = (int_expweight4) ? 
                        2.0*int_shell4->exponent(l)*ishl3expk : ishl3expk;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates(i,j,k,l, maxam1234);

          /* Generate the target integrals. */
          if (!maxam1234) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            intfunc brptr=build_routine[minam1][maxam12][minam3][maxam34][eAB];
            if (brptr == &BuildIntV3::impossible_integral) {
              ExEnv::errn() << scprintf("trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if (!(build.*brptr)()) {
              ExEnv::outn() << "build2e.cc: did not succeed in building all integrals"
                   << endl;
              abort();
              }
            }
          else {
            blockbuildprim(minam1,maxam12,minam3,maxam34);
            }

          /* Sum thru all possible contractions.
           * Throw out all unneeded contractions. */

  for (ci=0; ci<nc1; ci++) {
    int mlower = int_shell1->am(ci) + dam1;
    if (mlower < 0) continue;
    coef0 = int_shell1->coefficient_unnorm(ci,i)*c0scale;
    IntV3Arraydoublep2 ***e0f0_i = e0f0_con_ints_array[ci];
    for (cj=0; cj<nc2; cj++) {
      int mupper = mlower + int_shell2->am(cj) + dam2;
      if (mupper < mlower) continue;
      if (mlower < minam1) mlower = minam1;
      if (mupper > maxam12) mupper = maxam12;
      coef1 = int_shell2->coefficient_unnorm(cj,j)*coef0;
      IntV3Arraydoublep2 **e0f0_ij = e0f0_i[cj];
      for (ck=0; ck<nc3; ck++) {
        int nlower = int_shell3->am(ck) + dam3;
        if (nlower < 0) continue;
        coef2 = int_shell3->coefficient_unnorm(ck,k)*coef1;
        IntV3Arraydoublep2 *e0f0_ijk = e0f0_ij[ck];
        for (cl=0; cl<nc4; cl++) {
          int nupper = nlower + int_shell4->am(cl) + dam4;
          if (nupper < nlower) continue;
          if (nlower < minam3) nlower = minam3;
          if (nupper > maxam34) nupper = maxam34;
          coef3 = int_shell4->coefficient_unnorm(cl,l)*coef2;

          /* Contract the primitive target integrals. */
          if (i||j||k||l) {
            for (m=mlower; m<=mupper; m++) {
              int o;
              int sizec = contract_length(m,nlower,nupper);
              double *RESTRICT con_ints = e0f0_ijk[cl](m,nlower);
              bufferprim = build.int_v_list(m,nlower,0);
              /* Sum the integrals into the contracted integrals. */
#ifdef SUNMOS
              for (o=0; o < sizec; o++) {
                con_ints[o] += coef3 * bufferprim[o];
                }
#else
              for (o=sizec; o; o--) {
                *con_ints++ += coef3 * *bufferprim++;
                }
#endif
              }

            }
          else {
            for (m=mlower; m<=mupper; m++) {
              int o;
              int sizec = contract_length(m,nlower,nupper);
              double *RESTRICT con_ints = e0f0_ijk[cl](m,nlower);
              bufferprim = build.int_v_list(m,nlower,0);
              /* Write the integrals to the contracted integrals. */
#ifdef SUNMOS
              for (o=0; o < sizec; o++) {
                con_ints[o] = coef3 * bufferprim[o];
                }
#else
              for (o=sizec; o; o--) {
                *con_ints++ = coef3 * *bufferprim++;
                }
#endif
              }
            }
          }
        }
      }
    }


          }
        }
      }
    }
  }

/* This routine constructs intermediates needed for each quartet of
 * primitives.  It is given the total angular momentum as the argument
 * and requires that the global primitive offsets and other global
 * constants be initialized. */
void
Int2eV3::gen_prim_intermediates(int pr1, int pr2, int pr3, int pr4, int am)
{
  int i;
  double T;
  double pmq,pmq2;
  double AmB,AmB2;
  /* This is 2^(1/2) * pi^(5/4) */
  const double sqrt2pi54 = 5.9149671727956129;
  double conv_to_s;

  if (int_store2) {
    double *tmp;
    build.int_v_zeta12 = int_prim_zeta(opr1,opr2);
    build.int_v_zeta34 = int_prim_zeta(opr3,opr4);
    build.int_v_oo2zeta12 = int_prim_oo2zeta(opr1,opr2);
    build.int_v_oo2zeta34 = int_prim_oo2zeta(opr3,opr4);
    tmp = int_prim_p(opr1,opr2);
    build.int_v_p120 = *tmp++;
    build.int_v_p121 = *tmp++;
    build.int_v_p122 = *tmp;
    tmp = int_prim_p(opr3,opr4);
    build.int_v_p340 = *tmp++;
    build.int_v_p341 = *tmp++;
    build.int_v_p342 = *tmp;
    build.int_v_k12 = int_prim_k(opr1,opr2);
    build.int_v_k34 = int_prim_k(opr3,opr4);
    }
  else {
    build.int_v_zeta12 = int_shell1->exponent(pr1) + int_shell2->exponent(pr2);
    build.int_v_zeta34 = int_shell3->exponent(pr3) + int_shell4->exponent(pr4);
    build.int_v_oo2zeta12 = 1.0/build.int_v_zeta12;
    build.int_v_oo2zeta34 = 1.0/build.int_v_zeta34;
    build.int_v_p120 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r10
                         + int_shell2->exponent(pr2) * build.int_v_r20 );
    build.int_v_p121 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r11
                         + int_shell2->exponent(pr2) * build.int_v_r21 );
    build.int_v_p122 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r12
                         + int_shell2->exponent(pr2) * build.int_v_r22 );
    build.int_v_p340 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r30
                         + int_shell4->exponent(pr4) * build.int_v_r40 );
    build.int_v_p341 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r31
                         + int_shell4->exponent(pr4) * build.int_v_r41 );
    build.int_v_p342 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r32
                         + int_shell4->exponent(pr4) * build.int_v_r42 );

    /* Compute AmB^2 for shell 1 and 2. */
    AmB = build.int_v_r20 - build.int_v_r10;
    AmB2 = AmB*AmB;
    AmB = build.int_v_r21 - build.int_v_r11;
    AmB2 += AmB*AmB;
    AmB = build.int_v_r22 - build.int_v_r12;
    AmB2 += AmB*AmB;

    build.int_v_k12 =    sqrt2pi54
                 * build.int_v_oo2zeta12
                 * exp( -   int_shell1->exponent(pr1)*int_shell2->exponent(pr2)
                          * build.int_v_oo2zeta12
                          * AmB2 );

    /* Compute AmB^2 for shells 3 and 4. */
    AmB = build.int_v_r40 - build.int_v_r30;
    AmB2 = AmB*AmB;
    AmB = build.int_v_r41 - build.int_v_r31;
    AmB2 += AmB*AmB;
    AmB = build.int_v_r42 - build.int_v_r32;
    AmB2 += AmB*AmB;

    build.int_v_k34 =    sqrt2pi54
                 * build.int_v_oo2zeta34
                 * exp( -   int_shell3->exponent(pr3)*int_shell4->exponent(pr4)
                          * build.int_v_oo2zeta34
                          * AmB2 );

    build.int_v_oo2zeta12 *= 0.5;
    build.int_v_oo2zeta34 *= 0.5;
    }

  build.int_v_ooze = 1.0/(build.int_v_zeta12 + build.int_v_zeta34);

  build.int_v_W0 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p120
                         + build.int_v_zeta34 * build.int_v_p340 );
  build.int_v_W1 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p121
                         + build.int_v_zeta34 * build.int_v_p341 );
  build.int_v_W2 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p122
                         + build.int_v_zeta34 * build.int_v_p342 );

  pmq = build.int_v_p120 - build.int_v_p340;
  pmq2 = pmq*pmq;
  pmq = build.int_v_p121 - build.int_v_p341;
  pmq2 += pmq*pmq;
  pmq = build.int_v_p122 - build.int_v_p342;
  pmq2 += pmq*pmq;

  T =   build.int_v_zeta12
      * build.int_v_zeta34
      * build.int_v_ooze * pmq2;

  double *fjttable = fjt_->values(am,T);

  /* Convert the fjttable produced by int_fjt into the S integrals */
  conv_to_s = sqrt(build.int_v_ooze) * build.int_v_k12 * build.int_v_k34;
  for (i=0; i<=am; i++) {
    build.int_v_list(0,0,i)[0] =   fjttable[i] * conv_to_s;
    }

  }

/* This is like gen_prim_intermediates, except the normalization is
 * put into the ssss integrals. */
void
Int2eV3::gen_prim_intermediates_with_norm(int pr1, int pr2, int pr3, int pr4,
                                          int am, double norm)
{
  int i;
  double T;
  double pmq,pmq2;
  double AmB,AmB2;
  /* This is 2^(1/2) * pi^(5/4) */
  const double sqrt2pi54 = 5.9149671727956129;
  double conv_to_s;

  if (int_store2) {
    build.int_v_zeta12 = int_prim_zeta(opr1,opr2);
    build.int_v_zeta34 = int_prim_zeta(opr3,opr4);
    build.int_v_oo2zeta12 = int_prim_oo2zeta(opr1,opr2);
    build.int_v_oo2zeta34 = int_prim_oo2zeta(opr3,opr4);
    build.int_v_p120 = int_prim_p(opr1,opr2,0);
    build.int_v_p121 = int_prim_p(opr1,opr2,1);
    build.int_v_p122 = int_prim_p(opr1,opr2,2);
    build.int_v_p340 = int_prim_p(opr3,opr4,0);
    build.int_v_p341 = int_prim_p(opr3,opr4,1);
    build.int_v_p342 = int_prim_p(opr3,opr4,2);
    build.int_v_k12 = int_prim_k(opr1,opr2);
    build.int_v_k34 = int_prim_k(opr3,opr4);
    }
  else {
    build.int_v_zeta12 = int_shell1->exponent(pr1) + int_shell2->exponent(pr2);
    build.int_v_zeta34 = int_shell3->exponent(pr3) + int_shell4->exponent(pr4);
    build.int_v_oo2zeta12 = 1.0/build.int_v_zeta12;
    build.int_v_oo2zeta34 = 1.0/build.int_v_zeta34;
    build.int_v_p120 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r10
                         + int_shell2->exponent(pr2) * build.int_v_r20 );
    build.int_v_p121 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r11
                         + int_shell2->exponent(pr2) * build.int_v_r21 );
    build.int_v_p122 = build.int_v_oo2zeta12
                     * ( int_shell1->exponent(pr1) * build.int_v_r12
                         + int_shell2->exponent(pr2) * build.int_v_r22 );
    build.int_v_p340 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r30
                         + int_shell4->exponent(pr4) * build.int_v_r40 );
    build.int_v_p341 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r31
                         + int_shell4->exponent(pr4) * build.int_v_r41 );
    build.int_v_p342 = build.int_v_oo2zeta34
                     * ( int_shell3->exponent(pr3) * build.int_v_r32
                         + int_shell4->exponent(pr4) * build.int_v_r42 );

    /* Compute AmB^2 for shell 1 and 2. */
    AmB = build.int_v_r20 - build.int_v_r10;
    AmB2 = AmB*AmB;
    AmB = build.int_v_r21 - build.int_v_r11;
    AmB2 += AmB*AmB;
    AmB = build.int_v_r22 - build.int_v_r12;
    AmB2 += AmB*AmB;

    build.int_v_k12 =    sqrt2pi54
                 * build.int_v_oo2zeta12
                 * exp( -   int_shell1->exponent(pr1)*int_shell2->exponent(pr2)
                          * build.int_v_oo2zeta12
                          * AmB2 );

    /* Compute AmB^2 for shells 3 and 4. */
    AmB = build.int_v_r40 - build.int_v_r30;
    AmB2 = AmB*AmB;
    AmB = build.int_v_r41 - build.int_v_r31;
    AmB2 += AmB*AmB;
    AmB = build.int_v_r42 - build.int_v_r32;
    AmB2 += AmB*AmB;

    build.int_v_k34 =    sqrt2pi54
                 * build.int_v_oo2zeta34
                 * exp( -   int_shell3->exponent(pr3)*int_shell4->exponent(pr4)
                          * build.int_v_oo2zeta34
                          * AmB2 );

    build.int_v_oo2zeta12 *= 0.5;
    build.int_v_oo2zeta34 *= 0.5;
    }

  build.int_v_ooze = 1.0/(build.int_v_zeta12 + build.int_v_zeta34);

  build.int_v_W0 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p120
                         + build.int_v_zeta34 * build.int_v_p340 );
  build.int_v_W1 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p121
                         + build.int_v_zeta34 * build.int_v_p341 );
  build.int_v_W2 = build.int_v_ooze*(  build.int_v_zeta12 * build.int_v_p122
                         + build.int_v_zeta34 * build.int_v_p342 );

  pmq = build.int_v_p120 - build.int_v_p340;
  pmq2 = pmq*pmq;
  pmq = build.int_v_p121 - build.int_v_p341;
  pmq2 += pmq*pmq;
  pmq = build.int_v_p122 - build.int_v_p342;
  pmq2 += pmq*pmq;

  T =   build.int_v_zeta12
      * build.int_v_zeta34
      * build.int_v_ooze * pmq2;

  double *fjttable = fjt_->values(am,T);

  /* Convert the fjttable produced by int_fjt into the S integrals */
  conv_to_s = sqrt(build.int_v_ooze)
            * build.int_v_k12 * build.int_v_k34 * norm;
  for (i=0; i<=am; i++) {
    build.int_v_list(0,0,i)[0] =   fjttable[i] * conv_to_s;
    }

  }


/* This routine computes the shell intermediates. */
void
Int2eV3::gen_shell_intermediates(int sh1, int sh2, int sh3, int sh4)
{
  if (int_store1) {
    build.int_v_r10 = int_shell_r(osh1,0);
    build.int_v_r11 = int_shell_r(osh1,1);
    build.int_v_r12 = int_shell_r(osh1,2);
    build.int_v_r20 = int_shell_r(osh2,0);
    build.int_v_r21 = int_shell_r(osh2,1);
    build.int_v_r22 = int_shell_r(osh2,2);
    build.int_v_r30 = int_shell_r(osh3,0);
    build.int_v_r31 = int_shell_r(osh3,1);
    build.int_v_r32 = int_shell_r(osh3,2);
    build.int_v_r40 = int_shell_r(osh4,0);
    build.int_v_r41 = int_shell_r(osh4,1);
    build.int_v_r42 = int_shell_r(osh4,2);
    }
  else {
    build.int_v_r10 = pbs1_->r(pbs1_->shell_to_center(sh1),0);
    build.int_v_r11 = pbs1_->r(pbs1_->shell_to_center(sh1),1);
    build.int_v_r12 = pbs1_->r(pbs1_->shell_to_center(sh1),2);
    if (pbs2_ == 0) {
        build.int_v_r20 = 0.0;
        build.int_v_r21 = 0.0;
        build.int_v_r22 = 0.0;
      }
    else {
        build.int_v_r20 = pbs2_->r(pbs2_->shell_to_center(sh2),0);
        build.int_v_r21 = pbs2_->r(pbs2_->shell_to_center(sh2),1);
        build.int_v_r22 = pbs2_->r(pbs2_->shell_to_center(sh2),2);
      }
    build.int_v_r30 = pbs3_->r(pbs3_->shell_to_center(sh3),0);
    build.int_v_r31 = pbs3_->r(pbs3_->shell_to_center(sh3),1);
    build.int_v_r32 = pbs3_->r(pbs3_->shell_to_center(sh3),2);
    if (pbs4_ == 0) {
        build.int_v_r40 = 0.0;
        build.int_v_r41 = 0.0;
        build.int_v_r42 = 0.0;
      }
    else {
        build.int_v_r40 = pbs4_->r(pbs4_->shell_to_center(sh4),0);
        build.int_v_r41 = pbs4_->r(pbs4_->shell_to_center(sh4),1);
        build.int_v_r42 = pbs4_->r(pbs4_->shell_to_center(sh4),2);
      }
    }
  }

void
Int2eV3::blockbuildprim(int minam1,int maxam12,int minam3,int maxam34)
{
  int m, b;
  int l=maxam12+maxam34;

  // the (0,0,m) integrals have already been initialized

  // compute (0,b,m) integrals
  for (m=l-1; m>=0; m--) {
    int bmax = l-m;
    if (bmax>maxam34) bmax=maxam34;
    blockbuildprim_3(1,bmax,m);
    }

  // compute (a,b,m) integrals
  for (m=maxam12-1; m>=0; m--) {
    for (b=0; b<=maxam34; b++) {
// This is how the code was for a long while,
// but at some point it started giving the wrong
// answers and seems wrong from inspection.  Valgrind
// flags that uninitialized I10i integrals are being
// used, which results from amin > 1.  I have switched
// to the correctly behaving amin = 1.
//       int amin = minam1-m;
//       if (amin<1) amin=1;
//       int amax = maxam12-m;
//       blockbuildprim_1(amin,amax,b,m);
      int amax = maxam12-m;
      blockbuildprim_1(1,amax,b,m);
      }
    }
}

void
Int2eV3::blockbuildprim_1(int amin,int amax,int am34,int m)
{
  double *I00;
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a-1 0|c0](m) */
  double *I21; /* = [a-1 0|c0](m+1) */
  double *I31; /* = [a0|c-1 0](m+1) */
  int cartindex12;
  int cartindex34;
  int cartindex1234;
  int size34=0,size34m1;
  int i12, j12, k12;
  int i34, j34, k34;

  double ***vlist1;
  double **vlist10;
  double **vlist11;
  double ***vlist2;
  double **vlist20;

  vlist1 = build.int_v_list(amin-1);
  vlist10 = vlist1[am34];

  if (am34) {
    vlist11 = vlist1[am34-1];
    }

  if (amin>1) {
    vlist2 = build.int_v_list(amin-2);
    vlist20 = vlist2[am34];
    }

  /* The size of the am34 group of primitives. */
  size34 = INT_NCART_NN(am34);
  /* The size of the group of primitives with ang. mom. = am34 - 1 */
  size34m1 = INT_NCART_DEC(am34,size34);

  // Some local intermediates
  double half_ooze = 0.5 * build.int_v_ooze;
  double zeta34_ooze = build.int_v_zeta34 * build.int_v_ooze;
  double W0_m_p120 = build.int_v_W0 - build.int_v_p120;
  double p120_m_r10 = build.int_v_p120 - build.int_v_r10;
  double oo2zeta12 = build.int_v_oo2zeta12;
  double p121_m_r11 = build.int_v_p121 - build.int_v_r11;
  double W1_m_p121 = build.int_v_W1 - build.int_v_p121;
  double p122_m_r12 = build.int_v_p122 - build.int_v_r12;
  double W2_m_p122 = build.int_v_W2 - build.int_v_p122;

  stack_alignment_check(&half_ooze, "buildprim_1: half_ooze");

  for (int am12=amin; am12<=amax; am12++) {
    /* Construct the needed intermediate integrals. */
    double ***vlist0 = build.int_v_list(am12);
    double **vlist00 = vlist0[am34];
    I00 = vlist00[m];
    I10 = vlist10[m];
    I11 = vlist10[m+1];
    //I00 = build.int_v_list(am12,am34,m);
    //I10 = build.int_v_list(am12-1,am34,m);
    //I11 = build.int_v_list(am12-1,am34,m+1);
    if (am34) {
      I31 = vlist11[m+1];
      //I31 = build.int_v_list(am12 - 1, am34 - 1, m + 1);
      vlist11 = vlist0[am34-1];
      }
    if (am12>1) {
      I20 = vlist20[m];
      I21 = vlist20[m+1];
      //I20 = build.int_v_list(am12 - 2, am34, m);
      //I21 = build.int_v_list(am12 - 2, am34, m + 1);
      }
    vlist20 = vlist10;
    vlist10 = vlist00;

    /* Construct the new integrals. */
    cartindex12 = 0;
    cartindex1234 = 0;
    // the i12==0, k12==0, j12=am12 case (build on y)
    i12 = 0;
    j12 = am12;
    k12 = 0;

    int i12y1 = 0;  //= INT_CARTINDEX(am12-1,i12,j12-1);
    int i12y1s34 = i12y1*size34;
    int i12y1s34m1 = i12y1*size34m1;
    double *I10i = &I10[i12y1s34];
    double *I11i = &I11[i12y1s34];
    double *RESTRICT I00i = &I00[cartindex1234];
    if (j12==1) {
      for (cartindex34=0; cartindex34<size34; cartindex34++) {
        I00i[cartindex34]
          = I10i[cartindex34] * p121_m_r11
          + I11i[cartindex34] * W1_m_p121;
        }
      }
    else { // j12 > 1
      int i12y2s34 = 0; // = INT_CARTINDEX(am12-2,i12,j12-2)*size34;
      double *I20i = &I20[i12y2s34];
      double *I21i = &I21[i12y2s34];
      for (cartindex34=0; cartindex34<size34; cartindex34++) {
        I00i[cartindex34]
          = I10i[cartindex34] * p121_m_r11
          + I11i[cartindex34] * W1_m_p121
          + (j12 - 1) * oo2zeta12 * (I20i[cartindex34]
                                     - I21i[cartindex34] * zeta34_ooze);
        }
      }
    if (am34) {
      double *I31i = &I31[i12y1s34m1];
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        //note: k34 < am34-i34 instead of <= am34-i34, so j34 > 0
        int i34y1 = cartindex34-i34;//=INT_CARTINDEX(am34-1,i34,j34-1)
        j34 = am34 - i34;
        double j34_half_ooze = j34 * half_ooze;
        for (k34=0; k34<am34-i34; k34++) {

          I00i[cartindex34] +=  j34_half_ooze * I31i[i34y1];

          j34_half_ooze -= half_ooze;
          i34y1++;
          /* cartindex34 == INT_CARTINDEX(am34,i34,j34) */
          cartindex34++;
          }
        // increment cartindex34 here since the last k was skipped
        cartindex34++;
        }
      }
    cartindex12++;
    cartindex1234+=size34;

    // the i12==0, j12==am12-1, k12==1 case (build on z)
    i12 = 0;
    j12 = am12 - 1;
    k12 = 1;
    int i12z1 = 0;//= INT_CARTINDEX(am12-1,i12,j12);
    int i12z1s34 = i12z1*size34;
    int i12z1s34m1 = i12z1*size34m1;
    I10i = &I10[i12z1s34];
    I11i = &I11[i12z1s34];
    I00i = &I00[cartindex1234];
    for (cartindex34=0; cartindex34<size34; cartindex34++) {
      I00i[cartindex34]
        = I10i[cartindex34] * p122_m_r12
        + I11i[cartindex34] * W2_m_p122;
      }
    if (am34) {
      double *I31i = &I31[i12z1s34m1];
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        // skip k34 == 0
        cartindex34++;
        int i34z1 = cartindex34-i34-1;//=INT_CARTINDEX(am34-1,i34,j34)
        double k34_half_ooze = half_ooze;
        for (k34=1; k34<=am34-i34; k34++) {
          I00i[cartindex34] += k34_half_ooze * I31i[i34z1];
          k34_half_ooze += half_ooze;
          i34z1++;
          cartindex34++;
          }
        }
      }
    cartindex12++;
    cartindex1234+=size34;
    // the i12==0, j12==am12-k12, k12>1 case (build on z)
    double k12m1_oo2zeta12 = oo2zeta12;
    for (k12=2; k12<=am12-i12; k12++) {
      j12 = am12 - k12;
      i12z1 = cartindex12-i12-1;//=INT_CARTINDEX(am12-1,i12,j12);
      i12z1s34 = i12z1*size34;
      i12z1s34m1 = i12z1*size34m1;
      int i12z2s34 = (cartindex12-i12-i12-2)*size34;
      //=INT_CARTINDEX(am12-2,i12,j12)*size34;
      I10i = &I10[i12z1s34];
      I11i = &I11[i12z1s34];
      double *I20i = &I20[i12z2s34];
      double *I21i = &I21[i12z2s34];
      I00i = &I00[cartindex1234];
      for (cartindex34=0; cartindex34<size34; cartindex34++) {
        I00i[cartindex34]
          = I10i[cartindex34] * p122_m_r12
          + I11i[cartindex34] * W2_m_p122
          + k12m1_oo2zeta12 * (I20i[cartindex34]
                               - I21i[cartindex34] * zeta34_ooze);
        }
      if (am34) {
        double *I31i = &I31[i12z1s34m1];
        cartindex34 = 0;
        for (i34=0; i34<=am34; i34++) {
          // skip k34 == 0
          cartindex34++;
          int i34z1 = cartindex34-i34-1;//=INT_CARTINDEX(am34-1,i34,j34)
          double k34_half_ooze = half_ooze;
          for (k34=1; k34<=am34-i34; k34++) {
            I00i[cartindex34]
              +=  k34_half_ooze * I31i[i34z1];
            k34_half_ooze += half_ooze;
            i34z1++;
            cartindex34++;
            }
          }
        }
      cartindex12++;
      cartindex1234+=size34;
      k12m1_oo2zeta12 += oo2zeta12;
      }

    // the i12==1 case (build on x)
    i12 = 1;
    int i12x1 = cartindex12-am12-1;//=INT_CARTINDEX(am12-1,i12-1,am12-i12)
    int i12x1s34 = i12x1*size34;
    int i12x1s34m1 = i12x1*size34m1;
    I00i = &I00[cartindex1234];
    I10i = &I10[i12x1s34];
    I11i = &I11[i12x1s34];
    //for (k12=0; k12<=am12-i12; k12++)
    int k12_cartindex34;
    int nk12_size34 = am12*size34;
    for (k12_cartindex34=0; k12_cartindex34<nk12_size34; k12_cartindex34++) {
      *I00i++ = *I10i++ * p120_m_r10 + *I11i++ * W0_m_p120;
      }
    I00i = &I00[cartindex1234];
    if (am34) {
      double *I31i = &I31[i12x1s34m1];
      for (k12=0; k12<am12; k12++) {
        // skip over i34==0
        double *I00is=&I00i[am34+1];
        double i34_half_ooze = half_ooze;
        for (i34=1; i34<=am34; i34++) {
          for (k34=i34; k34<=am34; k34++) { // index_k34 = true_k34 + i34
            *I00is++ +=  i34_half_ooze * *I31i++;
            }
          i34_half_ooze += half_ooze;
          }
        I00i += size34;
        }
      }
    cartindex12 += am12;
    cartindex1234 += am12*size34;
    // the i12>1 case (build on x)
    if (am12<2) continue;
    double i12m1_oo2zeta12 = oo2zeta12;
    i12x1 = cartindex12-am12-1;
    i12x1s34 = i12x1*size34;
    i12x1s34m1 = i12x1*size34m1;
    int i12x2s34 = (cartindex12-am12-am12-1)*size34;
    I10i = &I10[i12x1s34];
    I11i = &I11[i12x1s34];
    double *I20i = &I20[i12x2s34];
    double *I21i = &I21[i12x2s34];
    I00i = &I00[cartindex1234];
    for (i12=2; i12<=am12; i12++) {
      int sizek12_size34 = (am12-i12+1)*size34;
      int k12_c34;
      for (k12_c34=0; k12_c34<sizek12_size34; k12_c34++) {
        *I00i++
          = *I10i++ * p120_m_r10
          + *I11i++ * W0_m_p120
          + i12m1_oo2zeta12 * (*I20i++
                               - *I21i++
                               * zeta34_ooze);
        }
      i12m1_oo2zeta12 += oo2zeta12;
      }
    if (am34) {
      double *I31i = &I31[i12x1s34m1];
      I00i = &I00[cartindex1234];
      for (i12=2; i12<=am12; i12++) {
        for (k12=0; k12<=am12-i12; k12++) {
          // skip over i34==0
          double *I00is=&I00i[am34+1];
          double i34_half_ooze = half_ooze;
          for (i34=1; i34<=am34; i34++) {
            for (k34=0; k34<=am34-i34; k34++) {
              *I00is++ += i34_half_ooze * *I31i++;
              }
            i34_half_ooze += half_ooze;
            }
          I00i += size34;
          }
        }
      }
    }
}

void
Int2eV3::blockbuildprim_3(int bmin,int bmax,int m)
{
  double *I00;
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a0|c-1 0](m) */
  double *I21; /* = [a0|c-1 0](m+1) */
  int ci34m1,ci34m2;
  int size34,size34m1,size34m2;
  int i34, k34;

  // These temporaries point to subblocks within the integrals arrays.
  double *I10o,*I11o,*I20o,*I21o;

  double ***vlist0;
  double **vlist01;
  double **vlist02;

  vlist0 = build.int_v_list(0);
  vlist01 = vlist0[bmin-1];
  if (bmin>1) {
    vlist02 = vlist0[bmin-2];
    }

  for (int am34=bmin; am34<=bmax; am34++) {

    /* Construct the needed intermediate integrals. */
    double **vlist00 = vlist0[am34];
    I00 = vlist00[m];
    I10 = vlist01[m];
    I11 = vlist01[m+1];
    //I00 = build.int_v_list(0, am34, m);
    //I10 = build.int_v_list(0, am34 - 1, m);
    //I11 = build.int_v_list(0, am34 - 1, m + 1);
    if (am34>1) {
      I20 = vlist02[m];
      I21 = vlist02[m+1];
      //I20 = build.int_v_list(0, am34 - 2, m);
      //I21 = build.int_v_list(0, am34 - 2, m + 1);
      }
    vlist02 = vlist01;
    vlist01 = vlist00;

    /* The size of the group of primitives with ang. mom. = am34 - 1 */
    size34 = INT_NCART_NN(am34);
    size34m1 = INT_NCART_DEC(am34,size34);
    size34m2 = INT_NCART(am34-2);

    // Useful constants
    double p340_m_r30 = build.int_v_p340 - build.int_v_r30;
    double W0_m_p340 = build.int_v_W0 - build.int_v_p340;
    double p341_m_r31 = build.int_v_p341 - build.int_v_r31;
    double W1_m_p341 = build.int_v_W1 - build.int_v_p341;
    double p342_m_r32 = build.int_v_p342 - build.int_v_r32;
    double W2_m_p342 = build.int_v_W2 - build.int_v_p342;
    double oo2zeta34 = build.int_v_oo2zeta34;
    double zeta12_ooze = build.int_v_zeta12 * build.int_v_ooze;

    stack_alignment_check(&p340_m_r30, "buildprim_3: p340_m_r30");

    /* Construct the new integrals. */
    double *RESTRICT I00o = I00; // points the current target integral
    I10o = I10;
    I11o = I11;
    //int cartindex34 = 0;
    // i34 == 0, k34 == 0, j34 = am34
    /* ------------------ Build from the y position. */
    /* I10 I11 and I21 */
    *I00o = *I10o * p341_m_r31 + *I11o * W1_m_p341;
    if (am34>1) {
      I20o = I20;
      I21o = I21;
      *I00o += (am34 - 1) * oo2zeta34 * (*I20o
                                         - *I21o * zeta12_ooze);
      }
    //cartindex34++;
    // i34 == 0, k34 >= 1
    // loop over a portion of the l=am34-1 integrals
    I00o = &I00o[1];
    for (ci34m1=0; ci34m1<am34; ci34m1++) {
      /* ------------------ Build from the z position. */
      //note: ci34m1 = cartindex34 - i34 - 1;//=INT_CARTINDEX(am34-1,i34,j34)
      /* I10 and I11 */
      I00o[ci34m1] = I10o[ci34m1] * p342_m_r32 + I11o[ci34m1] * W2_m_p342;
      }
    // skip over i34 == 0, k34 == 1
    //cartindex34++;
    // i34 == 0, k34 > 1
    I00o = &I00o[1];
    // loop over a portion of the l=am34-2 integrals
    double k34m1_oo2zeta34 = oo2zeta34;
    for (ci34m2=0; ci34m2<am34-1; ci34m2++) {
      //note: k34 = 2+ci34m2
      /* ------------------ Build from the z position. */
      /* I20 and I21 */
      I00o[ci34m2]
        +=  k34m1_oo2zeta34 * (I20o[ci34m2] - I21o[ci34m2] * zeta12_ooze);
      k34m1_oo2zeta34 += oo2zeta34;
      }
    //cartindex34+=am34-1;
    // i34 >= 1
    I00o = &I00o[am34-1];
    //note: ci34m1 = INT_CARTINDEX(am34-1,i34-1,j34)
    for (ci34m1=0; ci34m1<size34m1; ci34m1++) {
      /* I10 and I11 contrib */
      /* ------------------ Build from the x position. */
      I00o[ci34m1] = I10o[ci34m1] * p340_m_r30 + I11o[ci34m1] * W0_m_p340;
      }
    // skip past i34 == 1
    //cartindex34 += am34;
    // i34 > 1
    I00o = &I00o[am34];
    //note: ci34m2=INT_CARTINDEX(am34-2,i34-2,j34)
    ci34m2=0;
    double i34m1_oo2zeta34 = oo2zeta34;
    for (i34=2; i34<=am34; i34++) {
      for (k34=0; k34<=am34-i34; k34++) {
        /* I20 and I21 contrib */
        /* ------------------ Build from the x position. */
        I00o[ci34m2]
          +=  i34m1_oo2zeta34 * (I20o[ci34m2] - I21o[ci34m2] * zeta12_ooze);
        ci34m2++;
        }
      i34m1_oo2zeta34 += oo2zeta34;
      }
    //cartindex34 += size34m2;

    I00o = &I00o[size34m2];
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
