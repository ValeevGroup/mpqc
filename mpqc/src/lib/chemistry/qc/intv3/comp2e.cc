
#include <stdio.h>
#include <stdarg.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv3/tformv3.h>

static inline void
swtch(RefGaussianBasisSet &i,RefGaussianBasisSet &j)
{
  RefGaussianBasisSet tmp;
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
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
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
#ifdef TIMING
  char section[30];
#endif
  RefGaussianBasisSet pbs1=bs1_;
  RefGaussianBasisSet pbs2=bs2_;
  RefGaussianBasisSet pbs3=bs3_;
  RefGaussianBasisSet pbs4=bs4_;
  int size;
  int ii;
  int size1, size2, size3, size4;
  int tam1,tam2,tam3,tam4;
  int i,j,k,l;
  int ogc1,ogc2,ogc3,ogc4;
  int osh1,osh2,osh3,osh4;
  int sh1,sh2,sh3,sh4;
  int am1,am2,am3,am4,am12, am34;
  int minam1,minam2,minam3,minam4;
  int redundant_index;
  int e12,e13e24,e34;
  int i1,j1,k1;
  int i2,j2,k2;
  int i3,j3,k3;
  int i4,j4,k4;
  int p12,p34,p13p24;
  int eAB;

#if 0
  fprintf(stdout,"compute_erep: dam: (%d %d|%d %d)\n",
          dam1,dam2,dam3,dam4);
#endif

  /* Compute the offset shell numbers. */
  osh1 = *psh1 + bs1_shell_offset_;
  if (!int_unit2) osh2 = *psh2 + bs2_shell_offset_;
  osh3 = *psh3 + bs3_shell_offset_;
  if (!int_unit4) osh4 = *psh4 + bs4_shell_offset_;

  sh1 = *psh1;
  if (!int_unit2) sh2 = *psh2;
  sh3 = *psh3;
  if (!int_unit4) sh4 = *psh4;

  /* Test the arguments to make sure that they are sensible. */
  if (   sh1 < 0 || sh1 >= bs1_->nbasis()
      ||( !int_unit2 && (sh2 < 0 || sh2 >= bs2_->nbasis()))
      || sh3 < 0 || sh3 >= bs3_->nbasis()
      ||( !int_unit4 && (sh4 < 0 || sh4 >= bs4_->nbasis()))) {
    fprintf(stderr,"compute_erep has been incorrectly used\n");
    fprintf(stderr,"shells (bounds): %d (%d), %d (%d), %d (%d), %d (%d)\n",
            sh1,bs1_->nbasis()-1,
            sh2,bs2_->nbasis()-1,
            sh3,bs3_->nbasis()-1,
            sh4,bs4_->nbasis()-1);
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

#ifdef TIMING
  sprintf(section,"erep am=%d",am12+am34);
  tim_enter(section);
#endif

#if 0
  printf("on entry: (%d,%d,%d,%d) am=(%d,%d,%d,%d) perm = %d\n",
         *psh1, *psh2, *psh3, *psh4,
         am1,am2,am3,am4,
         !(INT_NOPERM&flags));
#endif

  /* Convert the integral to the most efficient form. */
  p12 = 0;
  p34 = 0;
  p13p24 = 0;

#undef OLD_PERMUTATION_ALGORITHM
#ifdef OLD_PERMUTATION_ALGORITHM
  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
    iswtch(&dam1,&dam2);
    iswtch(&minam1,&minam2);
    pswtch((void**)&int_shell1,(void**)&int_shell2);
    swtch(pbs1,pbs2);
    }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
    iswtch(&dam3,&dam4);
    iswtch(&minam3,&minam4);
    pswtch((void**)&int_shell3,(void**)&int_shell4);
    swtch(pbs3,pbs4);
    }
  if (!(int_unit2||int_unit4) && (osh1 == osh4) && (osh2 == osh3) && (osh1 != osh2)) {
    /* Don't make the permutation unless we won't override what was
     * decided above about p34. */
    if (am4 == am3) {
      p34 = 1;
      iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&dam3,&dam4);
      iswtch(&minam3,&minam4);
      pswtch((void**)&int_shell3,(void**)&int_shell4);
      swtch(pbs3,pbs4);
      }
    }
  if ((am3 > am1)||((am3 == am1)&&(am4 > am2))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&int_unit2,&int_unit4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    pswtch((void**)&int_shell1,(void**)&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell2,(void**)&int_shell4);
    swtch(pbs2,pbs4);
    }
  /* This tries to make centers A and B equivalent, if possible. */
  else if (  (am3 == am1)
           &&(am4 == am2)
           && !int_unit2
           && !int_unit4
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
    pswtch((void**)&int_shell1,(void**)&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell2,(void**)&int_shell4);
    swtch(pbs2,pbs4);
    }
#else /* OLD_PERMUTATION_ALGORITHM */
  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
    iswtch(&dam1,&dam2);
    iswtch(&minam1,&minam2);
    pswtch((void**)&int_shell1,(void**)&int_shell2);
    swtch(pbs1,pbs2);
    }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
    iswtch(&dam3,&dam4);
    iswtch(&minam3,&minam4);
    pswtch((void**)&int_shell3,(void**)&int_shell4);
    swtch(pbs3,pbs4);
    }
  if (!(int_unit2||int_unit4) && (osh1 == osh4) && (osh2 == osh3) && (osh1 != osh2)) {
    /* Don't make the permutation unless we won't override what was
     * decided above about p34. */
    if (am4 == am3) {
      p34 = 1;
      iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&dam3,&dam4);
      iswtch(&minam3,&minam4);
      pswtch((void**)&int_shell3,(void**)&int_shell4);
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
    pswtch((void**)&int_shell1,(void**)&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell2,(void**)&int_shell4);
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
    pswtch((void**)&int_shell1,(void**)&int_shell3);
    swtch(pbs1,pbs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell2,(void**)&int_shell4);
    swtch(pbs2,pbs4);
    }
#endif /* OLD_PERMUTATION_ALGORITHM */

  if (  int_unit2
        ||((pbs1 == pbs2)
          &&(pbs1->shell_to_center(sh1)==pbs2->shell_to_center(sh2)))) {
    eAB = 1;
    }
  else {
    eAB = 0;
    }

#if 0
  printf("perm info: p12 = %d, p34 = %d, p13p24 = %d, eAB = %d\n",
         p12,p34,p13p24,eAB);
#endif

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
  for (ii=size1=0; ii<int_shell1->ncontraction(); ii++)
    size1 += INT_NCART(int_shell1->am(ii)+dam1);
  for (ii=size2=0; ii<int_shell2->ncontraction(); ii++)
    size2 += INT_NCART(int_shell2->am(ii)+dam2);
  for (ii=size3=0; ii<int_shell3->ncontraction(); ii++)
    size3 += INT_NCART(int_shell3->am(ii)+dam3);
  for (ii=size4=0; ii<int_shell4->ncontraction(); ii++)
    size4 += INT_NCART(int_shell4->am(ii)+dam4);
  size = size1*size2*size3*size4;

  if (int_integral_storage) {
    if (dam1 || dam2 || dam3 || dam4) {
      fprintf(stderr,"cannot use integral storage and dam\n");
      fail();
      }
    if (    !int_unit2
         && !int_unit4
         && int_have_stored_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24))
      goto post_computation;
    }

  /* Buildam up on center 1 and 3. */
#if 0
  printf("C:");
#endif
  int_buildgcam(minam1,minam2,minam3,minam4,
                am1,am2,am3,am4,
                dam1,dam2,dam3,dam4,
                sh1,sh2,sh3,sh4, eAB);

  /* Begin loop over generalized contractions. */
  ogc1 = 0;
  for (i=0; i<int_shell1->ncontraction(); i++) {
    tam1 = int_shell1->am(i) + dam1;
    if (tam1 < 0) continue;
    ogc2 = 0;
    for (j=0; j<int_shell2->ncontraction(); j++) {
      tam2 = int_shell2->am(j) + dam2;
      if (tam2 < 0) continue;
      ogc3 = 0;
      for (k=0; k<int_shell3->ncontraction(); k++) {
        tam3 = int_shell3->am(k) + dam3;
        if (tam3 < 0) continue;
        ogc4 = 0;
        for (l=0; l<int_shell4->ncontraction(); l++) {
          tam4 = int_shell4->am(l) + dam4;
          if (tam4 < 0) continue;

  /* Shift angular momentum from 1 to 2 and from 3 to 4. */
  int_shiftgcam(i,j,k,l,tam1,tam2,tam3,tam4, eAB);

#if 0
  fprintf(stdout,"Top source int_con_ints_array[%d][%d][%d][%d]:\n",
          i,j,k,l);
  fprintf(stdout,"am = (%d,%d,%d,%d)  n = (%d,%d,%d,%d)\n",
         tam1+tam2,
         0,
         tam3+tam4,
         0,
         INT_NCART(tam1 + tam2),
         1,
         INT_NCART(tam3 + tam4),
         1
         );
  int_print_n(stdout,int_con_ints_array[i][j][k][l].dp
               [tam1 +tam2][0][tam3 +tam4][0],
               INT_NCART(tam1+tam2),
               1,
               INT_NCART(tam3+tam4),
               1,
               0,0,0);
  fprintf(stdout,"Targets:\n");
  fprintf(stdout,"am = (%d,%d,%d,%d)  n = (%d,%d,%d,%d)\n",
         tam1,
         tam2,
         tam3,
         tam4,
         INT_NCART(tam1),
         INT_NCART(tam2),
         INT_NCART(tam3),
         INT_NCART(tam4));
  int_print_n(stdout,int_con_ints_array[i][j][k][l].dp
                     [tam1]
                     [tam2]
                     [tam3]
                     [tam4],
               INT_NCART(tam1),
               INT_NCART(tam2),
               INT_NCART(tam3),
               INT_NCART(tam4),
               0,0,0);
#endif

  /* Place the integrals in the integral buffer. */
  /* If permute_ is not set, then repack the integrals while copying. */
  if ((!permute_)&&(p12||p34||p13p24)) {
    int newindex,pam1,pam2,pam3,pam4;
    int pi1,pj1,pk1;
    int pi2,pj2,pk2;
    int pi3,pj3,pk3;
    int pi4,pj4,pk4;
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
    FOR_CART(i1,j1,k1,tam1)
      FOR_CART(i2,j2,k2,tam2)
        FOR_CART(i3,j3,k3,tam3)
          FOR_CART(i4,j4,k4,tam4)
            pi1=i1; pj1=j1; pk1=k1;
            pi2=i2; pj2=j2; pk2=k2;
            pi3=i3; pj3=j3; pk3=k3;
            pi4=i4; pj4=j4; pk4=k4;
            if (p13p24) {
              iswtch(&pi1,&pi3);
              iswtch(&pj1,&pj3);
              iswtch(&pk1,&pk3);
              iswtch(&pi2,&pi4);
              iswtch(&pj2,&pj4);
              iswtch(&pk2,&pk4);
              }
            if (p34) {
              iswtch(&pi3,&pi4);
              iswtch(&pj3,&pj4);
              iswtch(&pk3,&pk4);
              }
            if (p12) {
              iswtch(&pi1,&pi2);
              iswtch(&pj1,&pj2);
              iswtch(&pk1,&pk2);
              }
            newindex =  (pogc1 + INT_CARTINDEX(pam1,pi1,pj1)) * psize234
                      + (pogc2 + INT_CARTINDEX(pam2,pi2,pj2)) * psize34
                      + (pogc3 + INT_CARTINDEX(pam3,pi3,pj3)) * psize4
                      + (pogc4 + INT_CARTINDEX(pam4,pi4,pj4));
#if 0
   printf("perm: int_buffer[%3d] += % f\n",newindex,
          int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index]);
#endif
            int_buffer[newindex]
              = int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index];
            redundant_index++;
            END_FOR_CART
          END_FOR_CART
        END_FOR_CART
      END_FOR_CART
    }
  else {
    int newindex;
    int size34 =  size3 * size4;
    int size234 = size2 * size34;
    double* redund_ints =
      int_con_ints_array[i][j][k][l].dp[tam1][tam2][tam3][tam4];
    redundant_index = 0;
#if 0
    printf("Using int_con_ints_array[%d][%d][%d][%d].dp[%d][%d][%d][%d]\n",
           i,j,k,l,tam1,tam2,tam3,tam4);
#endif
    FOR_CART(i1,j1,k1,tam1)
      FOR_CART(i2,j2,k2,tam2)
        FOR_CART(i3,j3,k3,tam3)
          FOR_CART(i4,j4,k4,tam4)
            newindex =  (ogc1 + INT_CARTINDEX(tam1,i1,j1)) * size234
                      + (ogc2 + INT_CARTINDEX(tam2,i2,j2)) * size34
                      + (ogc3 + INT_CARTINDEX(tam3,i3,j3)) * size4
                      + (ogc4 + INT_CARTINDEX(tam4,i4,j4));
#if 0
   printf("perm: int_buffer[%3d] += % f\n",newindex,
          int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index]);
#endif
            int_buffer[newindex] = redund_ints[redundant_index];
            redundant_index++;
            END_FOR_CART
          END_FOR_CART
        END_FOR_CART
      END_FOR_CART
    }

    /* End loop over generalized contractions. */
          ogc4 += INT_NCART(tam4);
          }
        ogc3 += INT_NCART(tam3);
        }
      ogc2 += INT_NCART(tam2);
      }
    ogc1 += INT_NCART(tam1);
    }
#if 0
  printf("\n");
#endif

  if (   !int_unit2
      && !int_unit4
      && int_integral_storage) {
      int_store_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24,size);
    }

  /* We branch here if an integral was precomputed and the int_buffer
   * has been already filled. */
  post_computation:

#if 0
  printf("before unpermute: am=(%d,%d,%d,%d)\n",am1,am2,am3,am4);
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
      pswtch((void**)&int_shell1,(void**)&int_shell3);
      swtch(pbs1,pbs3);
      pswtch((void**)&int_shell2,(void**)&int_shell4);
      swtch(pbs2,pbs4);
      iswtch(&int_expweight1,&int_expweight3);
      iswtch(&int_expweight2,&int_expweight4);
      }
    if (p34) {
      iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&am3,&am4);
      pswtch((void**)&int_shell3,(void**)&int_shell4);
      swtch(pbs3,pbs4);
      iswtch(&int_expweight3,&int_expweight4);
      }
    if (p12) {
      iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
      iswtch(&am1,&am2);
      pswtch((void**)&int_shell1,(void**)&int_shell2);
      swtch(pbs1,pbs2);
      iswtch(&int_expweight1,&int_expweight2);
      }
    }

  /* Transform to pure am (if requested in the centers structure). */
  if (!(flags&INT_NOPURE)) {
      intv3_transform_2e(int_buffer, int_buffer,
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
      fprintf(stderr,"nonredundant integrals cannot be generated\n");
      fail();
      }
    e12 = (int_unit2?0:(osh1 == osh2));
    e13e24 = ((osh1 == osh3)
              && ((int_unit2 && int_unit4)
                  || ((int_unit2||int_unit4)?0:(osh2 == osh4))));
    e34 = (int_unit4?0:(osh3 == osh4));
    nonredundant_erep(int_buffer,e12,e34,e13e24,
                           int_shell1->nfunction(),
                           int_shell2->nfunction(),
                           int_shell3->nfunction(),
                           int_shell4->nfunction(),
                           &redundant_offset,
                           &nonredundant_offset);
    }
    
#ifdef TIMING
  tim_exit(section);
#endif

#if 0
  printf("on exit: am=(%d,%d,%d,%d)\n",am1,am2,am3,am4);
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
  RefGaussianBasisSet cs[4];
  int sh[4];
  int n_unique;
  int i,j;
  GaussianShell *shell1,*shell2,*shell3,*shell4;
  RefGaussianBasisSet ucs[4]; /* The centers struct for the unique centers. */
  int ush[4];         /* The shells for the unique centers. */
  int unum[4];        /* The number of times that this unique center occurs. */
  int uam[4];         /* The total angular momentum on each unique center. */
  int am[4];
  int osh[4];
  int ncart;
  double *current_pure_buffer;

  cs[0] = bs1_;
  cs[1] = bs2_;
  cs[2] = bs3_;
  cs[3] = bs4_;

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

  /* This macro returns true if two shell centers are the same. */
#define SC(cs1,sh1,cs2,sh2) (((cs1)==(cs2))&&((cs1)->shell_to_center(sh1)==(cs1)->shell_to_center(sh2)))

  /* Build the list of unique centers structures and shells. */
  n_unique = 0;
  for (i=0; i<4; i++) {
    int unique = 1;
    for (j=0; j<n_unique; j++) {
      if (SC(ucs[j],ush[j],cs[i],sh[i])) {
        unique = 0;
        uam[j] += am[i];
        unum[j]++;
        break;
        }
      }
    if (unique) {
      ucs[n_unique] = cs[i];
      ush[n_unique] = sh[i];
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
    der_centers->num[der_centers->n] = ucs[i]->shell_to_center(ush[i]);
    der_centers->n++;

    for (j=0; j<4; j++) {
      if (SC(ucs[i],ush[i],cs[j],sh[j])) {
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
  der_centers->onum = ucs[omit]->shell_to_center(ush[omit]);

  /* Transform to pure am. */
  current_buffer = user_int_buffer;
  current_pure_buffer = user_int_buffer;
  for (i=0; i<3*der_centers->n; i++) {
      intv3_transform_2e(current_buffer, current_pure_buffer,
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
      fprintf(stderr,"nonredundant integrals cannot be generated (1der)\n");
      fail();
      }

    /* Shell equivalence information. */
    e12 = (osh[0] == osh[1]);
    e13e24 = ((osh[0] == osh[2]) && (osh[1] == osh[3]));
    e34 = (osh[2] == osh[3]);
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
  int oc1,oc2,oc3,oc4;
  int ii;
  int c1,c2,c3,c4;
  int i[4],j[4],k[4],am[4];
  int index;
  int sizem234,sizem34,sizem2,sizem3,sizem4;
  int sizep234,sizep34,sizep2,sizep3,sizep4;
  GaussianShell *shell1,*shell2,*shell3,*shell4;

  /* Set up pointers to the current shells. */
  shell1 = &bs1_->shell(*psh1);
  shell2 = &bs2_->shell(*psh2);
  shell3 = &bs3_->shell(*psh3);
  shell4 = &bs4_->shell(*psh4);

  if ((dercenter<0) || (dercenter > 3)) {
    fprintf(stderr,"illegal derivative center -- must be 0, 1, 2, or 3\n");
    fail();
    }

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
   /* The d/dx integrals */
  ALLDERLOOPS
    if (i[dercenter]>0) {
      buffer[index] -= i[dercenter] * int_buffer[
        (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0]-DCTEST(0),j[0])) * sizem234
       +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1]-DCTEST(1),j[1])) * sizem34
       +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2]-DCTEST(2),j[2])) * sizem4
       +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3]-DCTEST(3),j[3]))
       ];
      }
    index++;
    END_ALLDERLOOPS(-)

   /* The d/dy integrals */
  ALLDERLOOPS
    if (j[dercenter]>0) {
    buffer[index] -= j[dercenter] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0],j[0]-DCTEST(0))) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1],j[1]-DCTEST(1))) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2],j[2]-DCTEST(2))) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3],j[3]-DCTEST(3)))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

   /* The d/dz integrals */
  ALLDERLOOPS
    if (k[dercenter]>0) {
    buffer[index] -= k[dercenter] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0],j[0])) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1],j[1])) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2],j[2])) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3],j[3]))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

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
  /* The d/dx integrals */
  ALLDERLOOPS
#if 0
    tmp =int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
            ];
    if (1 || INT_NONZERO(tmp)) {
      printf("x: ((%d%d%d)(%d%d%d)(%d%d%d)(%d%d%d)) += ((%d%d%d)(%d%d%d)(%d%d%d)(%d%d%d)) (%f) (from %5d)\n",
             i[0],j[0],k[0],
             i[1],j[1],k[1],
             i[2],j[2],k[2],
             i[3],j[3],k[3],
             i[0]+DCTEST(0),j[0],k[0],
             i[1]+DCTEST(1),j[1],k[1],
             i[2]+DCTEST(2),j[2],k[2],
             i[3]+DCTEST(3),j[3],k[3],
             tmp,
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
             );
       }
#endif
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
            ];
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0],j[0]+DCTEST(0)))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1],j[1]+DCTEST(1)))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2],j[2]+DCTEST(2)))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3],j[3]+DCTEST(3)))
            ];
          index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCTEST(0),i[0],j[0])) * sizep234
              +(oc2 + INT_CARTINDEX(am[1]+DCTEST(1),i[1],j[1])) * sizep34
              +(oc3 + INT_CARTINDEX(am[2]+DCTEST(2),i[2],j[2])) * sizep4
              +(oc4 + INT_CARTINDEX(am[3]+DCTEST(3),i[3],j[3]))
              ];
          index++;
    END_ALLDERLOOPS(+)
  }

void
Int2eV3::nonredundant_erep(double *buffer, int e12, int e34, int e13e24,
                           int n1, int n2, int n3, int n4,
                           int *red_off, int *nonred_off)
{
  int redundant_index;
  int nonredundant_index;
  int i,j,k,l;

#ifdef I860
  if(mynode0() == -1 ) { printf("junk\n");}
#endif
  redundant_index = *red_off;
  nonredundant_index = *nonred_off;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      for (k=0; k<n3; k++) {
        for (l=0; l<n4; l++) {
          if (  (j<=INT_MAX2(e12,i,n2))
              &&(k<=INT_MAX3(e13e24,i,n3))
              &&(l<=INT_MAX4(e13e24,e34,i,j,k,n4))) {
            buffer[nonredundant_index] = buffer[redundant_index];
            nonredundant_index++;
            }
          redundant_index++;
          }
        }
      }
    }
  *red_off = redundant_index;
  *nonred_off = nonredundant_index;
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
      intv3_transform_2e(current_buffer, current_pure_buffer,
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
      fprintf(stderr,"nonredundant integrals cannot be generated (1der)\n");
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
    END_ALLDERLOOPS(+)

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
  END_ALLDERLOOPS(+)

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
  END_ALLDERLOOPS(+)

#if 0
  printf("after -DCT1 buffer[5] is %12.8f\n",buffer[5]);
#endif

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

#if 0
  printf("after +DCT1 buffer[5] is %12.8f\n",buffer[5]);
#endif

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

#if 0
  printf("after +DCT2 buffer[5] is %12.8f\n",buffer[5]);
#endif
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
