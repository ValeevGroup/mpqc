
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/fjt.h>
#include <chemistry/qc/intv3/utils.h>

#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/tformv3.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
}

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

  jmax1 = int_find_jmax(cs1);
  jmax2 = int_find_jmax(cs2);
  jmax = jmax1 + jmax2;

  fjt_ = new FJT(jmax + 2*order);

  nshell2 = int_find_ncartmax(cs1)*int_find_ncartmax(cs2);

  if (order == 0) {
    init_order = 0;
    scratchsize = nshell2;
    }
  else if (order == 1) {
    init_order = 1;
    scratchsize = nshell2*3;
    }
  else {
    fprintf(stderr,"int_initialize_1e: invalid order: %d\n",order);
    exit(1);
    }

#if 0
  printf("allocating %d doubles in init_1e\n",scratchsize);
#endif
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

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      cartesianbuffer[index] = comp_shell_overlap(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  intv3_transform_1e(cartesianbuffer, buff, shell1, shell2);
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
    fprintf(stderr,"int_shell_overlap: one electron routines are not init'ed\n");
    exit(1);
    }

  centers_t *dercs;
  if (idercs == 0) dercs = cs1;
  else dercs = cs2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

#if 0
  printf("zeroing %d*%d*3 elements of buff\n",ni,nj);
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
                                      centers_t *dercs, int centernum)
{
  accum_shell_1der(buff,ish,jsh,dercs,centernum,comp_shell_overlap);
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
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      exp1 = shell1->exp[i];
      exp2 = shell2->exp[j];
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
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_overlap(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      result += tmp;
      }
    }

  return result;
  }

/* Compute the overlap between two primitive functions. */
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
  oozeta = 1.0/(shell1->exp[prim1] + shell2->exp[prim2]);
  oo2zeta = 0.5*oozeta;
  AmB2 = 0.0;
  for (xyz=0; xyz<3; xyz++) {
    Pi = oozeta*(shell1->exp[prim1] * A[xyz] + shell2->exp[prim2] * B[xyz]);
    PmA[xyz] = Pi - A[xyz];
    PmB[xyz] = Pi - B[xyz];
    AmB = A[xyz] - B[xyz];
    AmB2 += AmB*AmB;
    }
  ss =   pow(3.141592653589793/(shell1->exp[prim1]+shell2->exp[prim2]),1.5)
       * exp(- oozeta * shell1->exp[prim1] * shell2->exp[prim2] * AmB2);
  return comp_prim_overlap(i1,j1,k1,i2,j2,k2);
  }

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

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  intv3_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

void
Int1eV3::int_accum_shell_kinetic(int ish, int jsh)
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int cart1,cart2;
  int index;
  int gc1,gc2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)
  intv3_accum_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the kinetic energy derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
void
Int1eV3::int_accum_shell_kinetic_1der(int ish, int jsh,
                                      centers_t *dercs, int centernum)
{
  accum_shell_1der(buff,ish,jsh,dercs,centernum,comp_shell_kinetic);
  }

/* This computes the basis function part of 
 * generic derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 * The function used to compute the nonderivative integrals is shell_function.
 */
void
Int1eV3::accum_shell_1der(double *buff, int ish, int jsh,
                          centers_t *dercs, int centernum,
                          double (Int1eV3::*shell_function)
                          (int,int,int,int,int,int,int,int))
{
  int i;
  int gc1,gc2;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index1,index2;
  double tmp[3];
  double *ctmp = cartesianbuffer;

  /* fprintf(stdout,"accum_shell_1der: working on ( %2d | %2d )",ish,jsh);
   * if (three_center) fprintf(stdout," third_centernum = %d\n",third_centernum);
   * else fprintf(stdout,"\n");
   */

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      if ((cs1==cs2)&&(c1==c2)) {
        if (    three_center
             && !((cs1==third_centers)&&(c1==third_centernum))
             && ((cs1==dercs)&&(c1==centernum))) {
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
      else if ((cs1==dercs)&&(c1==centernum)) {
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
      else if ((cs2==dercs)&&(c2==centernum)) {
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
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
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
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      oo2zeta_a = 0.5/shell1->exp[i];
      oo2zeta_b = 0.5/shell2->exp[j];
      oozeta = 1.0/(shell1->exp[i] + shell2->exp[j]);
      oo2zeta = 0.5*oozeta;
      xi = oozeta * shell1->exp[i] * shell2->exp[j];
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      /* The s integral kinetic energy. */
      ss =   pow(3.141592653589793/(shell1->exp[i]+shell2->exp[j]),1.5)
           * exp(- xi * AmB2);
      sTs =  ss
           * xi
           * (3.0 - 2.0 * xi * AmB2);
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_kinetic(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= shell1->exp[i];
      else if (exponent_weighted == 1) tmp *= shell2->exp[j];
      result += tmp;
      }
    }

  /* fprintf(stdout,"comp_shell_kinetic(%d,%d,%d,%d,%d,%d): result = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,result);
   */
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
                                      centers_t *dercs, int centernum)
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
                                          centers_t *dercs, int centernum)
{
  /* If both ish and jsh are not on the der center,
   * then there's no correction. */
  if (!(  (cs1==dercs)
        &&(cs2==dercs)
        &&(cs1->center_num[ish]==centernum)
        &&(cs2->center_num[jsh]==centernum))) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative for three equal
   * centers. */
  scale_shell_result = 1;
  result_scale_factor = -cs1->center[centernum].charge;
  C = cs1->center[centernum].r;
  accum_shell_efield(buff,ish,jsh);
  scale_shell_result = 0;

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the Hellman-Feynman part is computed.
 */
void
Int1eV3::int_accum_shell_nuclear_hf_1der(int ish, int jsh,
                                         centers_t *dercs, int centernum)
{

  /* If both ish and jsh are on the der center, then the contrib is zero. */
  if (  (cs1==dercs)
      &&(cs2==dercs)
      &&(cs1->center_num[ish]==centernum)
      &&(cs2->center_num[jsh]==centernum)) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative. */
  if (cs1 == dercs) {
    scale_shell_result = 1;
    result_scale_factor = -cs1->center[centernum].charge;
    C = cs1->center[centernum].r;
    accum_shell_efield(buff,ish,jsh);
    scale_shell_result = 0;
    }
  else if (cs2 == dercs) {
    scale_shell_result = 1;
    result_scale_factor = -cs2->center[centernum].charge;
    C = cs2->center[centernum].r;
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
                                            centers_t *dercs, int centernum)
{
  int i;

  /* Get the basis function part of the nuclear derivative. */
  three_center = 1;
  third_centers = cs1;
  for (i=0; i<cs1->n; i++) {
    third_centernum = i;
    C = cs1->center[i].r;
    scale_shell_result = 1;
    result_scale_factor = -cs1->center[i].charge;
    accum_shell_1der(buff,ish,jsh,dercs,centernum,comp_shell_nuclear);
    scale_shell_result = 0;
    }
  if (cs2!=cs1) {
    third_centers = cs2;
    for (i=0; i<cs2->n; i++) {
      third_centernum = i;
      C = cs2->center[i].r;
      scale_shell_result = 1;
      result_scale_factor = -cs2->center[i].charge;
      accum_shell_1der(buff,ish,jsh,dercs,centernum,comp_shell_nuclear);
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
  C = position;
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
    fprintf(stderr,"accum_shell_efield: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2);
      if (scale_shell_result) {
        for (i=0; i<3; i++) efield[i] *= result_scale_factor;
        }
      for (i=0; i<3; i++) tmp[i] = efield[i];
      tmp += 3;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the efield integrals between functions in two shells.
 * The result is placed in the buffer in the form bf1 x y z, bf2
 * x y z, etc.
 */
void
Int1eV3::efield(int ish, int jsh, double *position)
{
  scale_shell_result = 0;
  C = position;

  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    fprintf(stderr,"Int1eV3::efield one electron routines are not ready\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2);
      if (scale_shell_result) {
        for (i=0; i<3; i++) efield[i] *= result_scale_factor;
        }
      for (i=0; i<3; i++) tmp[i] = efield[i];
      tmp += 3;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  intv3_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
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
    fprintf(stderr,"int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      cartesianbuffer[index] = 0.0;
      /* Loop thru the centers on cs1. */
      for (i=0; i<cs1->n; i++) {
        C = cs1->center[i].r;
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * cs1->center[i].charge;
        }
      /* Loop thru the centers on cs2 if necessary. */
      if (cs2 != cs1) {
        for (i=0; i<cs2->n; i++) {
          C = cs2->center[i].r;
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                         * cs2->center[i].charge;
          }
        }
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  intv3_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
void
Int1eV3::int_accum_shell_point_charge(int ish, int jsh,
                                      int ncharge, double* charge,
                                      double** position)
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double tmp;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_pointcharge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      /* Loop thru the point charges. */
      tmp = 0.0;
      for (i=0; i<ncharge; i++) {
        C = position[i];
        tmp -=  comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * charge[i];
        }
      cartesianbuffer[index] = tmp;
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  intv3_accum_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
void
Int1eV3::point_charge(int ish, int jsh,
                                int ncharge, double* charge, double** position)
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;

  if (!(init_order >= 0)) {
    fprintf(stderr,"Int1eV3::point_charge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      cartesianbuffer[index] = 0.0;
      /* Loop thru the point charges. */
      for (i=0; i<ncharge; i++) {
        C = position[i];
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * charge[i];
        }
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  intv3_transform_1e(cartesianbuffer, buff, shell1, shell2);
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
    fprintf(stderr,"hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  index = 0;
  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      /* Loop thru the centers on cs1. */
      for (i=0; i<cs1->n; i++) {
        C = cs1->center[i].r;
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * cs1->center[i].charge;
        }
      /* Loop thru the centers on cs2 if necessary. */
      if (cs2 != cs1) {
        for (i=0; i<cs2->n; i++) {
          C = cs2->center[i].r;
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * cs2->center[i].charge;
          }
        }
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  intv3_transform_1e(cartesianbuffer, buff, shell1, shell2);
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
    fprintf(stderr,"int_shell_hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  centers_t *dercs;
  if (idercs == 0) dercs = cs1;
  else dercs = cs2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

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
    fprintf(stderr,"int_shell_kinetic: one electron routines are not init'ed\n");
    exit(1);
    }

  centers_t *dercs;
  if (idercs == 0) dercs = cs1;
  else dercs = cs2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

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
    fprintf(stderr,"int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  centers_t *dercs;
  if (idercs == 0) dercs = cs1;
  else dercs = cs2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

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
                                   centers_t *dercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear_hf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

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
                                      centers_t *dercs, int centernum)
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear_nonhf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

#if 0
  printf("int_shell_nuclear_nonhf_1der: zeroing %d doubles in buff\n",ni*nj*3);
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
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      zeta = shell1->exp[i] + shell2->exp[j];
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * 3.141592653589793/(shell1->exp[i]+shell2->exp[j])
           * exp(- oozeta * shell1->exp[i] * shell2->exp[j] * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am; k++) {
        fjttable_[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_nuclear(i1,j1,k1,i2,j2,k2,0);

      if (exponent_weighted == 0) tmp *= shell1->exp[i];
      else if (exponent_weighted == 1) tmp *= shell2->exp[j];

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

  /* fprintf(stdout,"  comp_prim_nuclear(%d,%d,%d,%d,%d,%d,%d) = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,m,result);
   */

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
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      zeta = shell1->exp[i] + shell2->exp[j];
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * 3.141592653589793/(shell1->exp[i]+shell2->exp[j])
           * exp(- oozeta * shell1->exp[i] * shell2->exp[j] * AmB2);

      /* The Fm(U) intermediates. */
      fjttable_ = fjt_->values(am+1,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am+1; k++) {
        fjttable_[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      for (xyz=0; xyz<3; xyz++) {
        result[xyz] +=  shell1->coef[gc1][i] * shell2->coef[gc2][j]
                      * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2,0);
        }
      }
    }

  for (xyz=0; xyz<3; xyz++) efield[xyz] = result[xyz];

  /* fprintf(stdout,"comp_shell_efield(%d,%d,%d,%d,%d,%d): % 12.8lf % 12.8lf % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,efield[0],efield[1],efield[2]);
   */
  }

double
Int1eV3::comp_prim_efield(int xyz, int i1, int j1, int k1,
                          int i2, int j2, int k2, int m)
{
  double result;

  /* if ((xyz != 0) || (i1 != 1)) return 0.0; */

  if (i1) {
    result  = PmA[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m);
    /* fprintf(stdout," % 12.8lf",result); */
    result -= PmC[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m+1);
    /* fprintf(stdout," % 12.8lf",result); */
    if (i1>1) result += oo2zeta * (i1-1)
                       * (  comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m+1));
    if (i2) result += oo2zeta * i2
                     * (  comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m+1));
    if (xyz==0) result += comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m+1);
    /* fprintf(stdout," % 12.8lf",result); */
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

  /* fprintf(stdout," % 12.8lf\n",result); */

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

  C = com;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2);
      for(mu=0; mu < 3; mu++) {
        cartesianbuffer[index] = dipole[mu];
        index++;
        }
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)
  intv3_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
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

  C = com;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2);
      for(mu=0; mu < 3; mu++) {
        cartesianbuffer[index] = dipole[mu];
        index++;
        }
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)
  intv3_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
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
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      exp1 = shell1->exp[i];
      exp2 = shell2->exp[j];
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
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j];
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
