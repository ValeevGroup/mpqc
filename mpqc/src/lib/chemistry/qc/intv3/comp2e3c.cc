
/* These routines compute two and three center electron repulsion
 * integrals. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/types.h>

#include <chemistry/qc/intv3/int2e.h>

void
Int2eV3::make_int_unit_shell()
{
  double *exp = new double[1];
  int *am = new int[1];
  int *pure = new int[1];
  double **c = new double*[1];
  *c = new double[1];

  exp[0] = 0.0;
  am[0] = 0;
  pure[0] = 0;
  c[0][0] = 1.0;

  int_unit_shell = new GaussianShell(1,1,exp,am,pure,c);
}

void
Int2eV3::delete_int_unit_shell()
{
  delete int_unit_shell;
  int_unit_shell = 0;
}

/* Compute a 2 center electron repulsion integral.  Electron 1 is in
 * shell psh1 and electron 2 is in psh2, that is (1 | 2).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
void
Int2eV3::erep_2center(int &psh1, int &psh2)
{
  RefGaussianBasisSet cs2 = bs2_;
  RefGaussianBasisSet cs4 = bs4_;
  int shd = 0x11111111; /* a dummy shell that will cause death if used */
  if (!int_unit_shell) make_int_unit_shell();
  bs2_ = 0;
  bs4_ = 0;
  int_unit2 = 1;
  int_unit4 = 1;
  erep(psh1,shd,psh2,shd);
  int_unit2 = 0;
  int_unit4 = 0;
  bs2_ = cs2;
  bs4_ = cs4;
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
      sizes[0] = bs1_->shell(shells[0]).nfunction();
      sizes[1] = bs3_->shell(shells[1]).nfunction();
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
      sizes[0] = bs1_->shell(shells[0]).nfunction();
      sizes[1] = bs3_->shell(shells[1]).nfunction();
      sizes[2] = bs4_->shell(shells[2]).nfunction();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
