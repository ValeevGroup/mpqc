
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <util/keyval/keyval.h>

#include <math/topology/point.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/cartiter.h>

#define MAX_NPRIM 20
#define MAX_NCON  10
#define MAX_AM    4

int
GaussianShell::values(const RefIntegral& ints,
                      const SCVector3& r, double* basis_values)
{
  return grad_values(ints, r, 0, basis_values);
}

// Returns a pointer to a vector of values of basis 
int
GaussianShell::grad_values(const RefIntegral& ints,
                           const SCVector3& r,
                           double* g_values,
                           double* basis_values) const
{

  // compute the maximum angular momentum component of the shell
  int maxam = max_am();
  if (g_values) maxam++;

  // check limitations
  if (nprim > MAX_NPRIM || ncon > MAX_NCON || maxam >= MAX_AM) {
      fprintf(stderr,"GaussianShell::grad_values: limit exceeded:\n");
      fprintf(stderr,"  ncon = %d (%d max) nprim = %d (%d max) maxam = %d (%d max)\n",
              ncon,MAX_NCON,nprim,MAX_NPRIM,maxam,MAX_AM-1);
      abort();
    }

  // loop variables
  int i,j;

  // precompute powers of x, y, and z
  double xs[MAX_AM];
  double ys[MAX_AM];
  double zs[MAX_AM];
  xs[0] = ys[0] = zs[0] = 1.0;
  if (maxam>0) {
      xs[1] = r[0];
      ys[1] = r[1];
      zs[1] = r[2];
    }
  for (i=2; i<=maxam; i++) {
      xs[i] = xs[i-1]*r[0];
      ys[i] = ys[i-1]*r[1];
      zs[i] = zs[i-1]*r[2];
    }

  // precompute r*r
  double r2;
  if (maxam<2) {
      r2 = 0.0;
      for (i=0; i<3; i++) {
          r2+=r[i]*r[i];
        }
    }
  else {
      r2 = xs[2] + ys[2] + zs[2];
    }

  // precompute exponentials
  double exps[MAX_NPRIM];
  for (i=0; i<nprim; i++) {
      exps[i]=::exp(-r2*exp[i]);
    }

  // precompute contractions over exponentials
  double precon[MAX_NCON];
  for (i=0; i<ncon; i++) {
      precon[i] = 0.0;
      for (j=0; j<nprim; j++) {
          precon[i] += coef[i][j] * exps[j];
        }
    }

  // precompute contractions over exponentials with exponent weighting
  double precon_g[MAX_NCON];
  if (g_values) {
      for (i=0; i<ncon; i++) {
          precon_g[i] = 0.0;
          for (j=0; j<nprim; j++) {
              precon_g[i] += exp[j] * coef[i][j] * exps[j];
            }
          precon_g[i] *= 2.0;
        }
    }

  // compute the shell values
  int i_basis=0;                // Basis function counter
  if (basis_values) {
      for (i=0; i<ncon; i++) {
          // handle s functions with a special case to speed things up
          if (l[i] == 0) {
              basis_values[i_basis] = precon[i];
              i_basis++;
            }
          else {
              CartesianIter *jp = ints->new_cartesian_iter(l[i]);
              CartesianIter& j = *jp;
              for (j.start(); j; j.next()) {
                  basis_values[i_basis] = xs[j.a()]*ys[j.b()]*zs[j.c()]
                                         *precon[i];
                  i_basis++;
                }
              delete jp;
            }
        }
    }

  // compute the gradient of the shell values
  if (g_values) {
      int i_grad=0;                // Basis function counter
      for (i=0; i<ncon; i++) {
          // handle s functions with a special case to speed things up
          if (l[i] == 0) {
              double norm_precon_g = precon_g[i];
              g_values[i_grad] = -xs[1]*norm_precon_g;
              i_grad++;
              g_values[i_grad] = -ys[1]*norm_precon_g;
              i_grad++;
              g_values[i_grad] = -zs[1]*norm_precon_g;
              i_grad++;
            }
          else {
              CartesianIter *jp = ints->new_cartesian_iter(l[i]);
              CartesianIter& j = *jp;
              for (j.start(); j; j.next()) {
                  double norm_precon = precon[i];
                  double norm_precon_g = precon_g[i];
                  g_values[i_grad] = - norm_precon_g
                    * xs[j.a()+1] * ys[j.b()] * zs[j.c()];
                  if (j.a()) g_values[i_grad] += j.a() * norm_precon
                    * xs[j.a()-1] * ys[j.b()] * zs[j.c()];
                  i_grad++;

                  g_values[i_grad] = - norm_precon_g
                    * xs[j.a()] * ys[j.b()+1] * zs[j.c()];
                  if (j.b()) g_values[i_grad] += j.b() * norm_precon
                    * xs[j.a()] * ys[j.b()-1] * zs[j.c()];
                  i_grad++;

                  g_values[i_grad] = - norm_precon_g
                    * xs[j.a()] * ys[j.b()] * zs[j.c()+1];
                  if (j.c()) g_values[i_grad] += j.c() * norm_precon
                    * xs[j.a()] * ys[j.b()] * zs[j.c()-1];
                  i_grad++;
                }
              delete jp;
            }
        }
    }

  return i_basis;
}
