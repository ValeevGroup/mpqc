
#include <stdio.h>
#include <stdlib.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsprnt.h>
}

/////////////////////////////////////////////////////////////////////////

/* Compute the number of shells, the shell offset for these centers,
 * and the arrays which map the number of the shell to the number
 * of the center and the number of the shell on that center. */
static int
shell_offset(centers_t *cs, int off)
{
  int ishell;
  int i,j;

  /* Set up the offset for the basis functions on this center. */
  cs->shell_offset = off;

  /* Compute the number of shells on this center. */
  cs->nshell = 0;
  for (i=0; i<cs->n; i++) {
    cs->nshell += cs->center[i].basis.n;
    }

  if (cs->center_num || cs->shell_num) {
    fprintf(stderr,"libint:shell_offset:center_num or shell_num is nonnull\n");
    print_centers(stderr,cs);
    exit(1);
    }

  /* Allocate storage for center_num and shell_num. */
  cs->center_num = (int *) malloc(sizeof(int)*cs->nshell);
  cs->shell_num = (int *) malloc(sizeof(int)*cs->nshell);
  if (!(cs->center_num && cs->shell_num)) {
    fprintf(stderr,"libint:shell_offset:problem allocating 3*%d integers\n",
            cs->nshell);
    exit(1);
    }

  /* Compute the center_num and shell_num arrays. */
  ishell = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      cs->center_num[ishell] = i;
      cs->shell_num[ishell] = j;
      ishell++;
      }
    }

  return off + cs->nshell;
  }

/* Compute function number array and the function offset
 * and the arrays which map the number of the shell to the number
 * of the first basis function in that that shell on that center. */
static int
func_offset(centers_t *cs, int off)
{
  int ishell;
  int i,j;

  /* Set up the offset for the basis functions on this center. */
  cs->func_offset = off;

  if (cs->func_num) {
    fprintf(stderr,"libint:func_offset:func_num is nonnull\n");
    print_centers(stderr,cs);
    exit(1);
    }

  /* Allocate storage for func_num. */
  cs->func_num = (int *) malloc(sizeof(int)*cs->nshell);
  if (!(cs->func_num)) {
    fprintf(stderr,"libint:shell_offset:problem allocating 3*%d integers\n",
            cs->nshell);
    exit(1);
    }

  if (!cs->nshell) return off;

  /* Compute the func_num array. */
  ishell = 0;
  cs->nfunc = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      shell_t *shell=&(cs->center[i].basis.shell[j]);
      cs->func_num[ishell] = cs->nfunc;
      cs->nfunc += shell->nfunc;
      ishell++;
      }
    }

  return off + cs->nfunc;
  }

/* Compute the primitive offsets. */
static int
prim_offset(centers_t *cs, int off)
{
  int i,j;

  /* Set up the offset for the basis functions on this center. */
  cs->prim_offset = off;

  /* Compute the number of primitives on this center. */
  cs->nprim = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      cs->nprim += cs->center[i].basis.shell[j].nprim;
      }
    }

  return off + cs->nprim;
  }

/* Compute the shell offset. */
static int
shell_offset(RefGaussianBasisSet cs, int off)
{
  return off + cs->nshell();
}

/* Compute the prim offset. */
static int
prim_offset(RefGaussianBasisSet cs, int off)
{
  return off + cs->nprimitive();
}

/* Compute the func offset. */
static int
func_offset(RefGaussianBasisSet cs, int off)
{
  return off + cs->nbasis();
}

/* Free storage for the offset arrays. */
static void
free_offsets(centers_t *cs)
{
  free(cs->center_num);
  free(cs->shell_num);
  free(cs->func_num);
  cs->center_num = NULL;
  cs->shell_num = NULL;
  cs->func_num = NULL;
  cs->shell_offset = 0;
  cs->prim_offset = 0;
  cs->func_offset = 0;
  cs->nshell = 0;
  cs->nprim = 0;
  cs->nfunc = 0;
  }

/////////////////////////////////////////////////////////////////////////

/* This initializes the offset arrays for one electron integrals. */
void
Int1eV3::int_initialize_offsets1(RefGaussianBasisSet cs1,
                                 RefGaussianBasisSet cs2)
{
  int shell_offset1;
  int prim_offset1;
  int func_offset1;

  /* Shell offset arrays. */
  bs1_shell_offset_ = 0;
  shell_offset1 = shell_offset(cs1,0);
  if (bs2_ != bs1_) {
      shell_offset(cs2,shell_offset1);
      bs2_shell_offset_ = shell_offset1;
    }
  else {
      bs2_shell_offset_ = bs1_shell_offset_;
    }

  /* Prim offset arrays. */
  bs1_prim_offset_ = 0;
  prim_offset1 = prim_offset(cs1,0);
  if (bs2_ != bs1_) {
      prim_offset(cs2,prim_offset1);
      bs2_prim_offset_ = prim_offset1;
    }
  else {
      bs2_prim_offset_ = bs1_prim_offset_;
    }

  /* Func offset arrays. */
  bs1_func_offset_ = 0;
  func_offset1 = func_offset(cs1,0);
  if (bs2_ != bs1_) {
      func_offset(cs2,func_offset1);
      bs2_func_offset_ = func_offset1;
    }
  else {
      bs2_func_offset_ = bs1_func_offset_;
    }
  }

/* This is called to free the offsets. */
void
Int1eV3::int_done_offsets1(RefGaussianBasisSet cs1, RefGaussianBasisSet cs2)
{
}

/* Initialize the offset arrays for two electron integrals. */
void
Int2eV3::int_initialize_offsets2(centers_t *cs1, centers_t *cs2,
                                 centers_t *cs3, centers_t *cs4)
{
  int shell_offset1;
  int shell_offset2;
  int shell_offset3;
  int shell_offset4;
  int prim_offset1;
  int prim_offset2;
  int prim_offset3;
  int prim_offset4;
  int func_offset1;
  int func_offset2;
  int func_offset3;
  int func_offset4;

  /* Shell offset arrays. */
  bs1_shell_offset_ = 0;

  shell_offset1 = shell_offset(cs1,0);
  if (bs2_ == bs1_) {
      shell_offset2 = shell_offset1;
      bs2_shell_offset_ = bs1_shell_offset_;
    }
  else {
      shell_offset2 = shell_offset(cs2,shell_offset1);
      bs2_shell_offset_ = shell_offset1;
    }

  if (bs3_ == bs1_) {
      shell_offset3 = shell_offset2;
      bs3_shell_offset_ = bs1_shell_offset_;
    }
  else if (bs3_ == bs2_) {
      shell_offset3 = shell_offset2;
      bs3_shell_offset_ = bs2_shell_offset_;
    }
  else {
      shell_offset3 = shell_offset(cs3,shell_offset2);
      bs3_shell_offset_ = shell_offset2;
    }

  if (bs4_ == bs1_) {
      shell_offset4 = shell_offset3;
      bs4_shell_offset_ = bs1_shell_offset_;
    }
  else if (bs4_ == bs2_) {
      shell_offset4 = shell_offset3;
      bs4_shell_offset_ = bs2_shell_offset_;
    }
  else if (bs4_ == bs3_) {
      shell_offset4 = shell_offset3;
      bs4_shell_offset_ = bs3_shell_offset_;
    }
  else {
      shell_offset4 = shell_offset(cs4,shell_offset3);
      bs4_shell_offset_ = shell_offset3;
    }

  /* Prim offset arrays. */
  bs1_prim_offset_ = 0;

  prim_offset1 = prim_offset(cs1,0);
  if (bs2_ == bs1_) {
      prim_offset2 = prim_offset1;
      bs2_prim_offset_ = bs1_prim_offset_;
    }
  else {
      prim_offset2 = prim_offset(cs2,prim_offset1);
      bs2_prim_offset_ = prim_offset1;
    }

  if (bs3_ == bs1_) {
      prim_offset3 = prim_offset2;
      bs3_prim_offset_ = bs1_prim_offset_;
    }
  else if (bs3_ == bs2_) {
      prim_offset3 = prim_offset2;
      bs3_prim_offset_ = bs2_prim_offset_;
    }
  else {
      prim_offset3 = prim_offset(cs3,prim_offset2);
      bs3_prim_offset_ = prim_offset2;
    }

  if (bs4_ == bs1_) {
      prim_offset4 = prim_offset3;
      bs4_prim_offset_ = bs1_prim_offset_;
    }
  else if (bs4_ == bs2_) {
      prim_offset4 = prim_offset3;
      bs4_prim_offset_ = bs2_prim_offset_;
    }
  else if (bs4_ == bs3_) {
      prim_offset4 = prim_offset3;
      bs4_prim_offset_ = bs3_prim_offset_;
    }
  else {
      prim_offset4 = prim_offset(cs4,prim_offset3);
      bs4_prim_offset_ = prim_offset3;
    }

  /* Func offset arrays. */
  bs1_func_offset_ = 0;

  func_offset1 = func_offset(cs1,0);
  if (bs2_ == bs1_) {
      func_offset2 = func_offset1;
      bs2_func_offset_ = bs1_func_offset_;
    }
  else {
      func_offset2 = func_offset(cs2,func_offset1);
      bs2_func_offset_ = func_offset1;
    }

  if (bs3_ == bs1_) {
      func_offset3 = func_offset2;
      bs3_func_offset_ = bs1_func_offset_;
    }
  else if (bs3_ == bs2_) {
      func_offset3 = func_offset2;
      bs3_func_offset_ = bs2_func_offset_;
    }
  else {
      func_offset3 = func_offset(cs3,func_offset2);
      bs3_func_offset_ = func_offset2;
    }

  if (bs4_ == bs1_) {
      func_offset4 = func_offset3;
      bs4_func_offset_ = bs1_func_offset_;
    }
  else if (bs4_ == bs2_) {
      func_offset4 = func_offset3;
      bs4_func_offset_ = bs2_func_offset_;
    }
  else if (bs4_ == bs3_) {
      func_offset4 = func_offset3;
      bs4_func_offset_ = bs3_func_offset_;
    }
  else {
      func_offset4 = func_offset(cs4,func_offset3);
      bs4_func_offset_ = func_offset3;
    }
  }

/* This is called to free the offsets. */
void
Int2eV3::int_done_offsets2(centers_t *cs1, centers_t *cs2,
                           centers_t *cs3, centers_t *cs4)
{
  free_offsets(cs1);
  if (!(cs2==cs1)) free_offsets(cs2);
  if (!((cs3==cs2)||(cs3==cs1))) free_offsets(cs3);
  if (!((cs4==cs3)||(cs4==cs2)||(cs4==cs1))) free_offsets(cs4);
  }
