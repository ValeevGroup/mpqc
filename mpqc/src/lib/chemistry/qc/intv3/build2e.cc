
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv3/macros.h>

#include <chemistry/qc/intv3/fjt.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/int2e.h>

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
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
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
  int e0f0_con_int_bufsize, con_int_bufsize;
  double *e0f0_con_int_buf, *con_int_buf;
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
  for (i=0; i<nc1; i++) {
    int tmp;
    jmax_for_con[i] = bs1_->max_am_for_contraction(i);
    if (  (bs2_ != bs1_)
        &&((tmp=bs2_->max_am_for_contraction(i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (bs3_ != bs1_) && (bs3_ != bs2_)
        &&((tmp=bs3_->max_am_for_contraction(i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (bs4_ != bs1_) && (bs4_ != bs2_) && (bs4_ != bs3_)
        &&((tmp=bs4_->max_am_for_contraction(i))>jmax_for_con[i]))
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
  if (  allocbn_int_array3(&inthave,"n1 n2 n3",am12+1,am34+1,am+1)
      ||allocbn_int_array3(&contract_length,"n1 n2 n3",am12+1,am34+1,am34+1)
      ||allocbn_doublep_array3(&build.int_v_list,"n1 n2 n3",am12+1,am34+1,am+1)) {
    fprintf(stderr,"problem allocating integral intermediates for");
    fprintf(stderr," am12 = %d, am34 = %d, and am = %d \n",am12,am34,am);
    fail();
    }

  /* Set all slots to NULL */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34; k++) {
        build.int_v_list.dp[i][j][k] = NULL;
        }
      }
    }

  for (i=0; i<=am12; i++) {
      for (j=0; j<=am34; j++) {
          for (k=0; k<=am34; k++) {
              contract_length.i[i][j][k] = 0;
              for (l=j; l<=k; l++) {
                  contract_length.i[i][j][k] += INT_NCART(i)*INT_NCART(l);
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

  int_v0_buf = (double*) malloc(sizeof(double)*int_v_bufsize);
  if (!int_v0_buf) {
    fprintf(stderr,"couldn't allocate all integral intermediates\n");
    fail();
    }
  add_store(int_v0_buf);
  int_v_buf = &int_v0_buf[int_v0_bufsize];

  /* Allocate storage for the needed slots. */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      build.int_v_list.dp[i][j][0] = int_v0_buf;
      int_v0_buf += INT_NCART(i)*INT_NCART(j);
      for (k=1; k<=am12+am34-i-j; k++) {
        build.int_v_list.dp[i][j][k] = int_v_buf;
        int_v_buf += INT_NCART(i)*INT_NCART(j);
        }
      }
    }


  /* Allocate storage for the contracted integrals (these are the output
   * of the build routines). */
  /* The ci, etc, indices refer to which of the pair of contraction
   * coefficients we are using. */
  e0f0_con_int_bufsize = 0;
  con_int_bufsize = 0;
  int_con_ints_array =
    (doublep_array4_t ****)malloc(sizeof(doublep_array4_t ***)*nc1);
  for (ci=0; ci<nc1; ci++) {
    int_con_ints_array[ci] =
      (doublep_array4_t ***)malloc(sizeof(doublep_array4_t **)*nc2);
    for (cj=0; cj<nc2; cj++) {
      int_con_ints_array[ci][cj] =
        (doublep_array4_t **)malloc(sizeof(doublep_array4_t *)*nc3);
      for (ck=0; ck<nc3; ck++) {
        int_con_ints_array[ci][cj][ck] =
          (doublep_array4_t *)malloc(sizeof(doublep_array4_t )*nc4);
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

  if (allocbn_doublep_array4(&int_con_ints_array[ci][cj][ck][cl],"n1 n2 n3 n4",
                             am12+1,am12+1,am34+1,am34+1)) {
    fprintf(stderr,"couldn't allocate contracted integral array\n");
    fail();
    }
  /* Count how much storage for the integrals is needed. */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am12; j++) {
      for (k=0; k<=am34; k++) {
        for (l=0; l<=am34; l++) {
          int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] = 0;
          }
        }
      }
    }
  for (i=0; i<=am12_for_con; i++) {
    for (j=0; j<=am12_for_con-i; j++) {
      for (k=0; k<=am34_for_con; k++) {
        for (l=0; l<=am34_for_con-k; l++) {
          if ((j==0)&&(l==0)) {
            int s =  INT_NCART(i)
                   * INT_NCART(j)
                   * INT_NCART(l)
                   * INT_NCART(k);
            e0f0_con_int_bufsize += s;
            con_int_bufsize += s;
             }
          else if ((ci==0)&&(cj==0)&&(ck==0)&&(cl==0)) {
            con_int_bufsize +=  INT_NCART(i)
                              * INT_NCART(j)
                              * INT_NCART(k)
                              * INT_NCART(l);
             }
          else if (int_con_ints_array[0][0][0][0].dp[i][j][k][l]) {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
              int_con_ints_array[0][0][0][0].dp[i][j][k][l];
             }
          else {
            con_int_bufsize +=  INT_NCART(i)
                              * INT_NCART(j)
                              * INT_NCART(k)
                              * INT_NCART(l);
             }
          }
        }
      }
    }
          }
        }
      }
    }
  e0f0_con_int_buf = (double*) malloc(con_int_bufsize*sizeof(double));
  used_storage_build_ += con_int_bufsize * sizeof(double);
  if (!e0f0_con_int_buf) {
    fprintf(stderr,"couldn't allocate contracted integral storage\n");
    fail();
    }
  add_store(e0f0_con_int_buf);
  con_int_buf = &e0f0_con_int_buf[e0f0_con_int_bufsize];
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

  for (i=0; i<=am12; i++) {
    for (j=0; j<=am12; j++) {
      for (k=0; k<=am34; k++) {
        for (l=0; l<=am34; l++) {
          int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] = 0;
          }
        }
      }
    }
  for (i=0; i<=am12_for_con; i++) {
    for (j=0; j<=am12_for_con-i; j++) {
      for (k=0; k<=am34_for_con; k++) {
        for (l=0; l<=am34_for_con-k; l++) {
/* If there are Pople style s=p shells and the shells are ordered
 * first s and then p and there are no p or d shells on the molecule,
 * then this algorithm would will allocate a little more storage
 * than needed.  General contraction should be ordered high to
 * low angular momentum for this reason. */
                /* Share storage for certain cases. */
          if ((j==0)&&(l==0)) {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]
                = e0f0_con_int_buf;
            e0f0_con_int_buf +=  INT_NCART(i)
                               * INT_NCART(j)
                               * INT_NCART(k)
                               * INT_NCART(l);
             }
           else if ((ci==0)&&(cj==0)&&(ck==0)&&(cl==0)) {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]
                = con_int_buf;
            con_int_buf +=  INT_NCART(i)
                          * INT_NCART(j)
                          * INT_NCART(k)
                          * INT_NCART(l);
               }
           else if (int_con_ints_array[0][0][0][0].dp[i][j][k][l]) {
             int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
               int_con_ints_array[0][0][0][0].dp[i][j][k][l];
             }
           else {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]
                = con_int_buf;
            con_int_buf +=  INT_NCART(i)
                          * INT_NCART(j)
                          * INT_NCART(k)
                          * INT_NCART(l);
             }
          }
        }
      }
    }
          }
        }
      }
    }

  /* This pointer is needed if int_buildam (without gc) is called. */
  int_con_ints = &int_con_ints_array[0][0][0][0];

  /* Initialize the build_routine array. */
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      for (k=0; k<4; k++) {
        for (l=0; l<4; l++) {
          for (m=0; m<2; m++) {
    build_routine[i][j][k][l][m] = BuildIntV3::impossible_integral;
            }
          }
        }
      }
    }

#define ASSIGN_BUILD(ii,j,k,l) \
  build_routine[ii][j][k][l][0]= BuildIntV3::i ## ii ## j ## k ## l ;\
  build_routine[ii][j][k][l][1]= BuildIntV3::i ## ii ## j ## k ## l ## eAB;
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
  }

void
Int2eV3::int_done_buildgc()
{
  int i,j,k;
  int ci,cj,ck,cl;

  used_storage_ -= used_storage_build_;
  used_storage_build_ = 0;

  free_int_array3(&inthave);
  free_int_array3(&contract_length);

  free_doublep_array3(&build.int_v_list);

  free_store();

  for (ci=0; ci<saved_ncon; ci++) {
    for (cj=0; cj<saved_ncon; cj++) {
      for (ck=0; ck<saved_ncon; ck++) {
        for (cl=0; cl<saved_ncon; cl++) {
          free_doublep_array4(&int_con_ints_array[ci][cj][ck][cl]);
          }
        free(int_con_ints_array[ci][cj][ck]);
        }
      free(int_con_ints_array[ci][cj]);
      }
    free(int_con_ints_array[ci]);
    }
  free(int_con_ints_array);

  }

/* add_store maintains a list of free storage allocated by int_init_buildgc */
void
Int2eV3::add_store(void *p)
{
  if (!store) {
    store = (store_list_t*) malloc(sizeof(store_list_t));
    assert(store);
    store->p = NULL;
    n_store_last = 0;
    }
  if (n_store_last >= STORAGE_CHUNK) {
    store_list_t* tmp = (store_list_t*) malloc(sizeof(store_list_t));
    assert(tmp);
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
  store = NULL;
  }

void
Int2eV3::_free_store(store_list_t* s, int n)
{
  int i;
  if (!s) return;
  for (i=0; i<n; i++) {
    free(s->data[i]);
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

#if 0
  fprintf(stdout,"min=(%d,%d,%d,%d) max=(%d,%d,%d,%d)\n",
          minam1,
          minam2,
          minam3,
          minam4,
          maxam1,
          maxam2,
          maxam3,
          maxam4);
#endif

  /* Compute the offset shell numbers. */
  osh1 = sh1 + bs1_shell_offset_;
  if (!int_unit2) osh2 = sh2 + bs2_shell_offset_;
  osh3 = sh3 + bs3_shell_offset_;
  if (!int_unit4) osh4 = sh4 + bs4_shell_offset_;

  nc1 = bs1_->shell(sh1).ncontraction();
  if (int_unit2) nc2 = 1;
  else nc2 = bs2_->shell(sh2).ncontraction();
  nc3 = bs3_->shell(sh3).ncontraction();
  if (int_unit4) nc4 = 1;
  else nc4 = bs4_->shell(sh4).ncontraction();

  /* Zero the target contracted integrals that the build routine
   * will accumulate into. */
  for (m=minam1; m<=maxam12; m++) {
    for (n=minam3; n<=maxam34; n++) {
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
      for (k=0; k<INT_NCART(m)*INT_NCART(n); k++) {
        int_con_ints_array[ci][cj][ck][cl].dp[m][0][n][0][k] = 0.0;
        }
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
  int have_all_ints;
  int i,j,k,l,m,n;
  int ci,cj,ck,cl;
  double *bufferprim;
  double *con_ints;

#if 0
  printf("not_gcs: %d%d%d%d\n",
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
    for (cj=0; cj<nc2; cj++) {
      int mupper = mlower + int_shell2->am(cj) + dam2;
      if (mupper < mlower) continue;
      if (mlower < minam1) mlower = minam1;
      if (mupper > maxam12) mupper = maxam12;
      for (ck=0; ck<nc3; ck++) {
        int nlower = int_shell3->am(ck) + dam3;
        if (nlower < 0) continue;
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
    if (int_store1) opr1 = int_shell_to_prim.i[osh1] + i;
    for (j=0; j<int_shell2->nprimitive(); j++) {
      double coef1;
      coef1 = int_shell2->coefficient_unnorm(cj,j);
      if (int_expweight2) coef1 *=  coef0
                                      * int_shell2->exponent(j);
      else                     coef1 *= coef0;
      if (int_store1) opr2 = int_shell_to_prim.i[osh2] + j;
      for (k=0; k<int_shell3->nprimitive(); k++) {
        double coef2;
        coef2 = int_shell3->coefficient_unnorm(ck,k);
        if (int_expweight3) coef2 *=  coef1
                                        * int_shell3->exponent(k);
        else                     coef2 *= coef1;
        if (int_store1) opr3 = int_shell_to_prim.i[osh3] + k;
        for (l=0; l<int_shell4->nprimitive(); l++) {
          double coef3;
          coef3 = int_shell4->coefficient_unnorm(cl,l);
          if (int_expweight4) coef3 *=  coef2
                                          * int_shell4->exponent(l);
          else                     coef3 *= coef2;
          if (int_store1) opr4 = int_shell_to_prim.i[osh4] + l;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates_with_norm(i,j,k,l, maxam12+maxam34,coef3);

          /* Generate the target integrals. */
          if ((maxam12 == 0) && (maxam34 == 0)) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            have_all_ints = 1;
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            if (build_routine[minam1]
                             [maxam12]
                             [minam3]
                             [maxam34][eAB]==BuildIntV3::impossible_integral){
              fprintf(stderr,"trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if ((build.*build_routine[minam1]
                                     [maxam12]
                                     [minam3]
                                     [maxam34][eAB])()) {
              /* Mark the integrals as computed. */
              have_all_ints = 1;
              }
            else {
              have_all_ints = 0;
              }
            }
          else if (MG > 0) {
            int backminam1  = MINA(minam1);
            int backmaxam12 = MINA(maxam12);
            int backminam3  = MINA(minam3);
            int backmaxam34 = MINA(maxam34);
            if ((backmaxam12==backmaxam34)&&(backminam1>backminam3))
              backminam3 = backminam1;
            /* We cannot build everything, so build what we can. */
            if (build_routine[backminam1]
                             [backmaxam12]
                             [backminam3]
                             [backmaxam34][eAB]
                == BuildIntV3::impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d\n",
                      backminam1,backmaxam12,backminam3,backmaxam34);
              }
            have_all_ints = 0;
            if ((build.*build_routine[backminam1]
                                     [backmaxam12]
                                     [backminam3]
                                     [backmaxam34][eAB])()) {

              /* Mark all the intermediates as being noncomputed. */
              init_inthave(maxam12,maxam34);

              /* Mark the integrals as computed. */
              for (m=backminam1; m<=backmaxam12; m++) {
                for (n=backminam3; n<=backmaxam34; n++) {
                  inthave.i[m][n][0] = 1;
                  }
                }
              }
            else {
              fprintf(stderr,"backup build routine not available\n");
              fail();
              }
            }
          else {
            init_inthave(maxam12,maxam34);
            have_all_ints = 0;
            }

          if (!have_all_ints) {
            for (m=minam1; m<=maxam12; m++) {
              if (m < int_shell1->am(ci)+dam1) continue;
              if (int_shell1->am(ci)+dam1+int_shell2->am(cj)+dam2
                  < m)
                continue;
              for (n=minam3; n<=maxam34; n++) {
                if (n < int_shell3->am(ck)+dam3) continue;
                if (int_shell3->am(ck)+dam3 +int_shell4->am(cl)+dam4
                    < n)
                  continue;

                buildprim(m, n, 0);
                }
              }
            have_all_ints = 1;
            }

          /* Contract the primitive target integrals. */
          /* Throw out all unneeded contractions. */
          for (m=mlower; m<=mupper; m++) {
            int o;
            int sizec = contract_length.i[m][nlower][nupper];
            con_ints = int_con_ints_array[ci][cj][ck][cl].dp[m][0][nlower][0];
            bufferprim = build.int_v_list.dp[m][nlower][0];

            for (o=sizec; o!=0; o--) {
              *con_ints++ += *bufferprim++;
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

#if 0
  print_int_array3(stdout,&inthave);
#endif

  }

void
Int2eV3::build_using_gcs(int nc1, int nc2, int nc3, int nc4,
                         int minam1, int minam3, int maxam12, int maxam34,
                         int dam1, int dam2, int dam3, int dam4, int eAB)
{
  int have_all_ints;
  int i,j,k,l,m,n;
  int ci,cj,ck,cl;
  int ist1,ist3;
  int nm3;
  int maxam1234=maxam12+maxam34;
  double coef0,coef1,coef2,coef3;
  double ishl1expi=1.0, ishl2expj=1.0, ishl3expk=1.0;
  double *bufferprim;
  double *con_ints;
  double c0scale;
  intfunc brptr=build_routine[minam1][maxam12][minam3][maxam34][eAB];

  /* Loop over the primitives. */
  for (i=0; i<int_shell1->nprimitive(); i++) {
    if (int_store1) opr1 = int_shell_to_prim.i[osh1] + i;
    if (int_expweight1) ishl1expi=2.0*int_shell1->exponent(i);

    for (j=0; j<int_shell2->nprimitive(); j++) {
      if (int_store1) opr2 = int_shell_to_prim.i[osh2] + j;
      ishl2expj = (int_expweight2) ? 
                        2.0*int_shell2->exponent(j)*ishl1expi : ishl1expi;

      for (k=0; k<int_shell3->nprimitive(); k++) {
        if (int_store1) opr3 = int_shell_to_prim.i[osh3] + k;
        ishl3expk = (int_expweight3) ? 
                        2.0*int_shell3->exponent(k)*ishl2expj : ishl2expj;

        for (l=0; l<int_shell4->nprimitive(); l++) {
          if (int_store1) opr4 = int_shell_to_prim.i[osh4] + l;
          c0scale = (int_expweight4) ? 
                        2.0*int_shell4->exponent(l)*ishl3expk : ishl3expk;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates(i,j,k,l, maxam1234);

          /* Generate the target integrals. */
          if (!maxam1234) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            have_all_ints = 1;
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            if (brptr == BuildIntV3::impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if ((build.*brptr)()) {
              /* Mark the integrals as computed. */
              have_all_ints = 1;
              }
            else {
              have_all_ints = 0;
              }
            }
          else if (MG > 0) {
            int backminam1  = MINA(minam1);
            int backmaxam12 = MINA(maxam12);
            int backminam3  = MINA(minam3);
            int backmaxam34 = MINA(maxam34);
            intfunc brptr2;

            if ((backmaxam12==backmaxam34)&&(backminam1>backminam3))
              backminam3 = backminam1;
            brptr2=build_routine
                       [backminam1][backmaxam12][backminam3][backmaxam34][eAB];

            /* We cannot build everything, so build what we can. */
            if (brptr2 == BuildIntV3::impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d\n",
                      backminam1,backmaxam12,backminam3,backmaxam34);
              }
            have_all_ints = 0;
            if ((build.*brptr2)()) {
              /* Mark all the intermediates as being noncomputed. */
              init_inthave(maxam12,maxam34);

              /* Mark the integrals as computed. */
              for (m=backminam1; m<=backmaxam12; m++) {
                for (n=backminam3; n<=backmaxam34; n++) {
                  inthave.i[m][n][0] = 1;
                  }
                }
              }
            else {
              fprintf(stderr,"backup build routine not available\n");
              fail();
              }
            }
          else {
            init_inthave(maxam12,maxam34);
            have_all_ints = 0;
            }

          /* Sum thru all possible contractions.
           * Throw out all unneeded contractions. */

          if (!have_all_ints) {
            for (m=minam1; m<=maxam12; m++) {
              for (n=minam3; n<=maxam34; n++) {
                bufferprim = buildprim(m, n, 0);
                }
              }
            have_all_ints = 1;
            }

  for (ci=0; ci<nc1; ci++) {
    int mlower = int_shell1->am(ci) + dam1;
    if (mlower < 0) continue;
    coef0 = int_shell1->coefficient_unnorm(ci,i)*c0scale;
    for (cj=0; cj<nc2; cj++) {
      int mupper = mlower + int_shell2->am(cj) + dam2;
      if (mupper < mlower) continue;
      if (mlower < minam1) mlower = minam1;
      if (mupper > maxam12) mupper = maxam12;
      coef1 = int_shell2->coefficient_unnorm(cj,j)*coef0;
      for (ck=0; ck<nc3; ck++) {
        int nlower = int_shell3->am(ck) + dam3;
        if (nlower < 0) continue;
        coef2 = int_shell3->coefficient_unnorm(ck,k)*coef1;
        for (cl=0; cl<nc4; cl++) {
          int nupper = nlower + int_shell4->am(cl) + dam4;
          if (nupper < nlower) continue;
          if (nlower < minam3) nlower = minam3;
          if (nupper > maxam34) nupper = maxam34;
          coef3 = int_shell4->coefficient_unnorm(cl,l)*coef2;

          /* Contract the primitive target integrals. */
          for (m=mlower; m<=mupper; m++) {
            int o;
            int sizec = contract_length.i[m][nlower][nupper];
            con_ints = int_con_ints_array[ci][cj][ck][cl].dp[m][0][nlower][0];
            bufferprim = build.int_v_list.dp[m][nlower][0];

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

  if (int_store2 && !int_unit2 && !int_unit4) {
    double *tmp;
    build.int_v_zeta12 = int_prim_zeta.d[opr1][opr2];
    build.int_v_zeta34 = int_prim_zeta.d[opr3][opr4];
    build.int_v_oo2zeta12 = int_prim_oo2zeta.d[opr1][opr2];
    build.int_v_oo2zeta34 = int_prim_oo2zeta.d[opr3][opr4];
    tmp = int_prim_p.d[opr1][opr2];
    build.int_v_p120 = *tmp++;
    build.int_v_p121 = *tmp++;
    build.int_v_p122 = *tmp;
    tmp = int_prim_p.d[opr3][opr4];
    build.int_v_p340 = *tmp++;
    build.int_v_p341 = *tmp++;
    build.int_v_p342 = *tmp;
    build.int_v_k12 = int_prim_k.d[opr1][opr2];
    build.int_v_k34 = int_prim_k.d[opr3][opr4];
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
    build.int_v_list.dp[0][0][i][0] =   fjttable[i] * conv_to_s;
#if 0
    fprintf(stdout,"build.int_v_list.dp[0][0][%d][0] = %lf\n",
            i,build.int_v_list.dp[0][0][i][0]);
#endif
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
    build.int_v_zeta12 = int_prim_zeta.d[opr1][opr2];
    build.int_v_zeta34 = int_prim_zeta.d[opr3][opr4];
    build.int_v_oo2zeta12 = int_prim_oo2zeta.d[opr1][opr2];
    build.int_v_oo2zeta34 = int_prim_oo2zeta.d[opr3][opr4];
    build.int_v_p120 = int_prim_p.d[opr1][opr2][0];
    build.int_v_p121 = int_prim_p.d[opr1][opr2][1];
    build.int_v_p122 = int_prim_p.d[opr1][opr2][2];
    build.int_v_p340 = int_prim_p.d[opr3][opr4][0];
    build.int_v_p341 = int_prim_p.d[opr3][opr4][1];
    build.int_v_p342 = int_prim_p.d[opr3][opr4][2];
    build.int_v_k12 = int_prim_k.d[opr1][opr2];
    build.int_v_k34 = int_prim_k.d[opr3][opr4];
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
    build.int_v_list.dp[0][0][i][0] =   fjttable[i] * conv_to_s;
#if 0
    fprintf(stdout,"build.int_v_list.dp[0][0][%d][0] = %lf\n",
            i,build.int_v_list.dp[0][0][i][0]);
#endif
    }

  }


/* This routine computes the shell intermediates. */
void
Int2eV3::gen_shell_intermediates(int sh1, int sh2, int sh3, int sh4)
{
  if (int_store1 && !int_unit2 && !int_unit4) {
    build.int_v_r10 = int_shell_r.d[osh1][0];
    build.int_v_r11 = int_shell_r.d[osh1][1];
    build.int_v_r12 = int_shell_r.d[osh1][2];
    build.int_v_r20 = int_shell_r.d[osh2][0];
    build.int_v_r21 = int_shell_r.d[osh2][1];
    build.int_v_r22 = int_shell_r.d[osh2][2];
    build.int_v_r30 = int_shell_r.d[osh3][0];
    build.int_v_r31 = int_shell_r.d[osh3][1];
    build.int_v_r32 = int_shell_r.d[osh3][2];
    build.int_v_r40 = int_shell_r.d[osh4][0];
    build.int_v_r41 = int_shell_r.d[osh4][1];
    build.int_v_r42 = int_shell_r.d[osh4][2];
    }
  else {
    build.int_v_r10 = bs1_->r(bs1_->shell_to_center(sh1),0);
    build.int_v_r11 = bs1_->r(bs1_->shell_to_center(sh1),1);
    build.int_v_r12 = bs1_->r(bs1_->shell_to_center(sh1),2);
    if (int_unit2) {
        build.int_v_r20 = build.int_v_r10;
        build.int_v_r21 = build.int_v_r11;
        build.int_v_r22 = build.int_v_r12;
      }
    else {
        build.int_v_r20 = bs2_->r(bs2_->shell_to_center(sh2),0);
        build.int_v_r21 = bs2_->r(bs2_->shell_to_center(sh2),1);
        build.int_v_r22 = bs2_->r(bs2_->shell_to_center(sh2),2);
      }
    build.int_v_r30 = bs3_->r(bs3_->shell_to_center(sh3),0);
    build.int_v_r31 = bs3_->r(bs3_->shell_to_center(sh3),1);
    build.int_v_r32 = bs3_->r(bs3_->shell_to_center(sh3),2);
    if (int_unit4) {
        build.int_v_r40 = build.int_v_r30;
        build.int_v_r41 = build.int_v_r31;
        build.int_v_r42 = build.int_v_r32;
      }
    else {
        build.int_v_r40 = bs4_->r(bs4_->shell_to_center(sh4),0);
        build.int_v_r41 = bs4_->r(bs4_->shell_to_center(sh4),1);
        build.int_v_r42 = bs4_->r(bs4_->shell_to_center(sh4),2);
      }
    }
  }

/* This builds up the primitive integrals of the type [x0|y0](m). */
double *
Int2eV3::buildprim(int am12, int am34, int m)
{
  double *buffer;

  /* Is this no integral? */
  if ((am12 < 0) || (am34 < 0)) return NULL;

  /* Is this integral on the list of computed integrals? */
  if (inthave.i[am12][am34][m]) return build.int_v_list.dp[am12][am34][m];

  /* Find the preallocated storage for the integrals. */
  buffer = build.int_v_list.dp[am12][am34][m];

  /* Should we build on center 1 or center 3. */
  if (choose_center(am12,am34,m) == 1) {
    /* Build on 1. */
    buildprim_1(buffer,am12,am34,m);
    }
  else {
    /* Build on 3. */
    buildprim_3(buffer,am12,am34,m);
    }

  /* Put the integrals in the list of computed integrals. */
  inthave.i[am12][am34][m] = 1;

#if 0
  fprintf(stdout,"buildprim: (%d 0|%d 0)(%d):\n",am12,am34,m);
  int_print_n(stdout,buffer
              ,INT_NCART(am12)
              ,1
              ,INT_NCART(am34)
              ,1
              ,0,0,0);
#endif

  return buffer;
  }

/* I00 will be made [a+1 0|b 0](m) */
void
Int2eV3::buildprim_1(double *I00, int am12, int am34, int m)
{
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a-1 0|c0](m) */
  double *I21; /* = [a-1 0|c0](m+1) */
  double *I31; /* = [a0|c-1 0](m+1) */
  int cartindex12;
  int cartindex34;
  int cartindex1234;
  int size34,size34m1;
  int i12, j12, k12;
  int i34, j34, k34;

  /* Construct the needed intermediate integrals. */
  I10 = buildprim(am12 - 1, am34, m);
  I11 = buildprim(am12 - 1, am34, m + 1);
  I20 = buildprim(am12 - 2, am34, m);
  I21 = buildprim(am12 - 2, am34, m + 1);
  I31 = buildprim(am12 - 1, am34 - 1, m + 1);

  /* The size of the am34 group of primitives. */
  size34 = INT_NCART(am34);
  /* The size of the group of primitives with ang. mom. = am34 - 1 */
  size34m1 = INT_NCART(am34-1);

  /* Construct the new integrals. */
  cartindex12 = 0;
  cartindex1234 = 0;
  for (i12=0; i12<=am12; i12++) {
    for (k12=0; k12<=am12-i12; k12++) {
      j12 = am12 - i12 - k12;
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        for (k34=0; k34<=am34-i34; k34++) {
          j34 = am34 - i34 - k34;

          /* I10 I11 I20 I21 and I31 contrib */

          /* ------------------ Build from the x position. */
          if (i12) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
                * ( build.int_v_p120 - build.int_v_r10)
               + I11[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
                 * ( build.int_v_W0 - build.int_v_p120);
            if (i12 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (i12 - 1) * build.int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
                      * build.int_v_zeta34 * build.int_v_ooze);
              }
            if (i34) {
              /* I31 */
              I00[cartindex1234]
               +=  i34 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34-1,j34)];
              }
            }
          /* ------------------ Build from the y position. */
          else if (j12) {
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12,j12-1)*size34 + cartindex34]
                * ( build.int_v_p121 - build.int_v_r11)
               + I11[INT_CARTINDEX(am12-1,i12,j12-1)*size34 + cartindex34]
                 * ( build.int_v_W1 - build.int_v_p121);
            if (j12 > 1) {
              I00[cartindex1234]
               +=  (j12 - 1) * build.int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12,j12-2)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12,j12-2)*size34 + cartindex34]
                      * build.int_v_zeta34 * build.int_v_ooze);
              }
            if (j34) {
              I00[cartindex1234]
               +=  j34 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12-1)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34-1)];
              }
            }
          /* ------------------ Build from the z position. */
          else if (k12) {
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12,j12)*size34 + cartindex34]
                * ( build.int_v_p122 - build.int_v_r12)
               + I11[INT_CARTINDEX(am12-1,i12,j12)*size34 + cartindex34]
                 * ( build.int_v_W2 - build.int_v_p122);
            if (k12 > 1) {
              I00[cartindex1234]
               +=  (k12 - 1) * build.int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12,j12)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12,j12)*size34 + cartindex34]
                      * build.int_v_zeta34 * build.int_v_ooze);
              }
            if (k34) {
              I00[cartindex1234]
               +=  k34 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34)];
              }
            }

          /* cartindex34 == INT_CARTINDEX(am34,i34,j34) */
          cartindex34++;
          cartindex1234++;
          }
        }
      /* cartindex12 == INT_CARTINDEX(am12,i12,j12) */
      cartindex12++;
      }
    }

  }


/* I00 will be made [a 0|b+1 0](m) */
void
Int2eV3::buildprim_3(double *I00, int am12, int am34, int m)
{
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a0|c-1 0](m) */
  double *I21; /* = [a0|c-1 0](m+1) */
  double *I31; /* = [a-1 0|c0](m+1) */
  int cartindex12;
  int cartindex34;
  int cartindex1234;
  int size34m1,size34m2;
  int i12, j12, k12;
  int i34, j34, k34;

  /* Construct the needed intermediate integrals. */
  I10 = buildprim(am12, am34 - 1, m);
  I11 = buildprim(am12, am34 - 1, m + 1);
  I20 = buildprim(am12, am34 - 2, m);
  I21 = buildprim(am12, am34 - 2, m + 1);
  I31 = buildprim(am12 - 1, am34 - 1, m + 1);

  /* The size of the group of primitives with ang. mom. = am34 - 1 */
  size34m1 = INT_NCART(am34-1);
  size34m2 = INT_NCART(am34-2);

  /* Construct the new integrals. */
  cartindex12 = 0;
  cartindex1234 = 0;
  for (i12=0; i12<=am12; i12++) {
    for (k12=0; k12<=am12-i12; k12++) {
      j12 = am12 - i12 - k12;
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        for (k34=0; k34<=am34-i34; k34++) {
          j34 = am34 - i34 - k34;

          /* I10 I11 I20 I21 and I31 contrib */

          /* ------------------ Build from the x position. */
          if (i34) {
            /* I10 and I11 */
            I00[cartindex1234]
        /*
         *    = I10[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
         *      * ( build.int_v_p120 - build.int_v_r10)
         *     + I11[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
         *       * ( build.int_v_W0 - build.int_v_p120);
         */
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34-1,j34)]
                * ( build.int_v_p340 - build.int_v_r30)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34-1,j34)]
                 * ( build.int_v_W0 - build.int_v_p340);
            if (i34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
         /*
          *    +=  (i12 - 1) * build.int_v_oo2zeta12
          *      * (  I20[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
          *         - I21[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
          *           * build.int_v_zeta34 * build.int_v_ooze);
          */
               +=  (i34 - 1) * build.int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34-2,j34)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34-2,j34)]
                      * build.int_v_zeta12 * build.int_v_ooze);
              }
            if (i12) {
              /* I31 */
              I00[cartindex1234]
        /*  
         *     +=  (i34 - 1) * 0.5 * build.int_v_ooze
         *        * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
         *              + INT_CARTINDEX(am34-1,i34-1,j34)];
         */
               +=  i12 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34-1,j34)];
              }
            }
          /* ------------------ Build from the y position. */
          else if (j34) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34-1)]
                * ( build.int_v_p341 - build.int_v_r31)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34-1)]
                 * ( build.int_v_W1 - build.int_v_p341);
            if (j34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (j34 - 1) * build.int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34-2)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34-2)]
                      * build.int_v_zeta12 * build.int_v_ooze);
              }
            if (j12) {
              /* I31 */
              I00[cartindex1234]
               +=  j12 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12-1)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34-1)];
              }
            }
          /* ------------------ Build from the z position. */
          else if (k34) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34)]
                * ( build.int_v_p342 - build.int_v_r32)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34)]
                 * ( build.int_v_W2 - build.int_v_p342);
            if (k34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (k34 - 1) * build.int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34)]
                      * build.int_v_zeta12 * build.int_v_ooze);
              }
            if (k12) {
              /* I31 */
              I00[cartindex1234]
               +=  k12 * 0.5 * build.int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34)];
              }
            }

          /* cartindex34 == INT_CARTINDEX(am34,i34,j34) */
          cartindex34++;
          cartindex1234++;
          }
        }
      /* cartindex12 == INT_CARTINDEX(am12,i12,j12) */
      cartindex12++;
      }
    }

  }

/* Initialize the list of integrals which have been precomputed
 * to "not computed" (=0). */
void
Int2eV3::init_inthave(int am12, int am34)
{
  int i,j,k;

#if 0 /* This cleans up inthave if the whole thing is to be printed */
  zero_int_array3(&inthave);
#endif

  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34-i-j; k++) {
        if ((i==0)&&(j==0)) inthave.i[i][j][k] = 1;
        else inthave.i[i][j][k] = 0;
        }
      }
    }
  }


/* Determines which center it is best to build upon.
 * (This can be improved, perhaps.) */
int
Int2eV3::choose_center(int am12, int am34, int m)
{
  int need1 = 0;
  int need3 = 0;

  if (am12==0) return 3;
  if (am34==0) return 1;

  if (!inthave.i[am12-1][am34][m]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34);
    }
  if (!inthave.i[am12-1][am34][m+1]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34);
    }
  if ((am12>1) && (!inthave.i[am12-2][am34][m])) {
    need1 += INT_NCART(am12-2)*INT_NCART(am34);
    }
  if ((am12>1) && (!inthave.i[am12-2][am34][m+1])) {
    need1 += INT_NCART(am12-2)*INT_NCART(am34);
    }
  if (!inthave.i[am12-1][am34-1][m+1]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34-1);
    }

  if (!inthave.i[am12][am34-1][m]) {
    need3 += INT_NCART(am12)*INT_NCART(am34-1);
    }
  if (!inthave.i[am12][am34-1][m+1]) {
    need3 += INT_NCART(am12)*INT_NCART(am34-1);
    }
  if ((am34>1) && (!inthave.i[am12][am34-2][m])) {
    need3 += INT_NCART(am12)*INT_NCART(am34-2);
    }
  if ((am34>1) && (!inthave.i[am12][am34-2][m+1])) {
    need3 += INT_NCART(am12)*INT_NCART(am34-2);
    }
  if (!inthave.i[am12-1][am34-1][m+1]) {
    need3 += INT_NCART(am12-1)*INT_NCART(am34-1);
    }

  if (need1 <= need3) return 1;
  return 3;
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
