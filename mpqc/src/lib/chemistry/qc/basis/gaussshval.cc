//
// gaussshval.cc
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
#include <util/keyval/keyval.h>

#include <math/topology/point.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>

#define MAX_NPRIM 20
#define MAX_NCON  10
#define MAX_AM    6

int
GaussianShell::values(CartesianIter **civec, SphericalTransformIter **sivec,
                      const SCVector3& r, double* basis_values)
{
  return grad_values(civec, sivec, r, 0, basis_values);
}

// Returns a pointer to a vector of values of basis 
int
GaussianShell::grad_values(CartesianIter **civec,
                           SphericalTransformIter **sivec,
                           const SCVector3& r,
                           double* g_values,
                           double* basis_values) const
{

  // compute the maximum angular momentum component of the shell
  int maxam = max_am();
  if (g_values) maxam++;

  // check limitations
  if (nprim > MAX_NPRIM || ncon > MAX_NCON || maxam >= MAX_AM) {
      cerr << node0 << indent
           << "GaussianShell::grad_values: limit exceeded:\n"
           << indent
           << scprintf(
               "ncon = %d (%d max) nprim = %d (%d max) maxam = %d (%d max)\n",
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
          else if (!puream[i]) {
              CartesianIter *jp = civec[l[i]];
              CartesianIter& j = *jp;
              for (j.start(); j; j.next()) {
                  basis_values[i_basis] = xs[j.a()]*ys[j.b()]*zs[j.c()]
                                         *precon[i];
                  i_basis++;
                }
            }
          else {
              double cart_basis_values[((MAX_AM+1)*(MAX_AM+2))/2];
              CartesianIter *jp = civec[l[i]];
              CartesianIter& j = *jp;
              int i_cart = 0;
              for (j.start(); j; j.next()) {
                  cart_basis_values[i_cart] = xs[j.a()]*ys[j.b()]*zs[j.c()]
                                             *precon[i];
                  i_cart++;
                }
              SphericalTransformIter *ti = sivec[l[i]];
              int n = ti->n();
              memset(&basis_values[i_basis], 0, sizeof(double)*n);
              for (ti->start(); ti->ready(); ti->next()) {
                  basis_values[i_basis + ti->pureindex()]
                      += ti->coef() * cart_basis_values[ti->cartindex()];
                }
              i_basis += n;
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
          else if (!puream[i]) {
              CartesianIter *jp = civec[l[i]];
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
            }
          else {
              double cart_g_values[3*((MAX_AM+1)*(MAX_AM+2))/2];
              CartesianIter *jp = civec[l[i]];
              CartesianIter& j = *jp;
              int i_cart = 0;
              for (j.start(); j; j.next()) {
                  double norm_precon = precon[i];
                  double norm_precon_g = precon_g[i];
                  cart_g_values[i_cart] = - norm_precon_g
                    * xs[j.a()+1] * ys[j.b()] * zs[j.c()];
                  if (j.a()) cart_g_values[i_cart] += j.a() * norm_precon
                    * xs[j.a()-1] * ys[j.b()] * zs[j.c()];
                  i_cart++;

                  cart_g_values[i_cart] = - norm_precon_g
                    * xs[j.a()] * ys[j.b()+1] * zs[j.c()];
                  if (j.b()) cart_g_values[i_cart] += j.b() * norm_precon
                    * xs[j.a()] * ys[j.b()-1] * zs[j.c()];
                  i_cart++;

                  cart_g_values[i_cart] = - norm_precon_g
                    * xs[j.a()] * ys[j.b()] * zs[j.c()+1];
                  if (j.c()) cart_g_values[i_cart] += j.c() * norm_precon
                    * xs[j.a()] * ys[j.b()] * zs[j.c()-1];
                  i_cart++;
                }
              SphericalTransformIter *ti = sivec[l[i]];
              int n = ti->n();
              memset(&g_values[i_grad], 0, sizeof(double)*n*3);
              for (ti->start(); ti->ready(); ti->next()) {
                  double coef = ti->coef();
                  int pi = ti->pureindex();
                  int ci = ti->cartindex();
                  for (int xyz=0; xyz<3; xyz++) {
                      g_values[i_grad + pi*3 + xyz]
                          += coef * cart_g_values[ci*3 + xyz];
                    }
                }
              i_grad += 3*n;
            }
        }
    }

  return i_basis;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
