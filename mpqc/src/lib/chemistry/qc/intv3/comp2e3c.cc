
/* These routines compute two and three center electron repulsion
 * integrals. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/types.h>

#include <chemistry/qc/intv3/int2e.h>

void
Int2eV3::make_int_unit_shell()
{
  int_unit_shell = (shell_t*) malloc(sizeof(shell_t));
  int_unit_shell->nprim = 1;
  int_unit_shell->ncon = 1;
  int_unit_shell->nfunc = 1;
  int_unit_shell->exp = (double*) malloc(sizeof(double));
  int_unit_shell->exp[0] = 0.0;
  int_unit_shell->type = (shell_type_t*) malloc(sizeof(shell_type_t));
  int_unit_shell->type[0].am = 0;
  int_unit_shell->type[0].puream = 0;
  int_unit_shell->coef = (double**) malloc(sizeof(double*));
  int_unit_shell->coef[0] = (double*) malloc(sizeof(double));
  int_unit_shell->coef[0][0] = 1.0;
}

void
Int2eV3::delete_int_unit_shell()
{
  if (!int_unit_shell) return;
  free(int_unit_shell->coef[0]);
  free(int_unit_shell->coef);
  free(int_unit_shell->type);
  free(int_unit_shell->exp);
  free(int_unit_shell);
  int_unit_shell = 0;
}

/* Compute a 2 center electron repulsion integral.  Electron 1 is in
 * shell psh1 and electron 2 is in psh2, that is (1 | 2).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
void
Int2eV3::erep_2center(int &psh1, int &psh2)
{
  centers_t *cs2 = int_cs2;
  centers_t *cs4 = int_cs4;
  int shd = 0x11111111; /* a dummy shell that will cause death if used */
  if (!int_unit_shell) make_int_unit_shell();
  int_cs2 = (centers_t*) 0x0;
  int_cs4 = (centers_t*) 0x0;
  int_unit2 = 1;
  int_unit4 = 1;
  erep(psh1,shd,psh2,shd);
  int_unit2 = 0;
  int_unit4 = 0;
  int_cs2 = cs2;
  int_cs4 = cs4;
}

/* This is an alternate interface to int_erep2.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
void
Int2eV3::erep_2center(int *shells, int  *sizes)
{
  erep_2center(shells[0],shells[1]);
  if (sizes) {
      sizes[0] = INT_SH(int_cs1,shells[0]).nfunc;
      sizes[1] = INT_SH(int_cs3,shells[1]).nfunc;
    }
}


/* Computes a 3 center two electron integral.  Electron 1 is in psh1
 * and electron 2 is in psh2 and psh3, that is (1 | 2 3).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
void
Int2eV3::erep_3center(int &psh1, int &psh2, int &psh3)
{
  int shd = 0x11111111; /* a dummy shell that will cause death if used */
  int oldperm = permute();
  if (!int_unit_shell) make_int_unit_shell();
  int_unit2 = 1;
  erep(psh1,shd,psh2,psh3);
  int_unit2 = 0;
  set_permute(oldperm);
}

/* This is an alternate interface to int_erep3.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
void
Int2eV3::erep_3center(int *shells, int  *sizes)
{
  erep_3center(shells[0],shells[1],shells[2]);
  if (sizes) {
      sizes[0] = INT_SH(int_cs1,shells[0]).nfunc;
      sizes[1] = INT_SH(int_cs3,shells[1]).nfunc;
      sizes[2] = INT_SH(int_cs4,shells[2]).nfunc;
    }
}
