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

#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>

using namespace std;
using namespace sc;

#define MAX_NPRIM 20
#define MAX_NCON  10
#define MAX_AM    8

int
GaussianShell::values(CartesianIter **civec, SphericalTransformIter **sivec,
                      const SCVector3& r, double* basis_values)
{
  return hessian_values(civec, sivec, r, 0, 0, basis_values);
}

int
GaussianShell::grad_values(CartesianIter **civec,
                           SphericalTransformIter **sivec,
                           const SCVector3& r,
                           double* g_values,
                           double* basis_values) const
{
  return hessian_values(civec, sivec, r, 0, g_values, basis_values);
}

int
GaussianShell::hessian_values(CartesianIter **civec,
                           SphericalTransformIter **sivec,
                           const SCVector3& r,
                           double* h_values,
                           double* g_values,
                           double* basis_values) const
{

  // compute the maximum angular momentum component of the shell
  int maxam = max_am();
  if (g_values || h_values) maxam++;
  if (h_values) maxam++;

  // check limitations
  if (nprimitive() > MAX_NPRIM || ncontraction() > MAX_NCON || maxam >= MAX_AM) {
      std::ostringstream oss;
      oss
           << "GaussianShell::grad_values: limit exceeded:\n"
           << indent
           << scprintf(
               "ncon = %d (%d max) nprim = %d (%d max) maxam = %d (%d max)\n",
               ncontraction(),MAX_NCON,nprimitive(),MAX_NPRIM,maxam,MAX_AM-1);
      ExEnv::out0() << oss.str();
      throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
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
  for (i=0; i<nprimitive(); i++) {
      exps[i]=::exp(-r2*exp[i]);
    }

  // precompute contractions over exponentials
  double precon[MAX_NCON];
  for (i=0; i<ncontraction(); i++) {
      precon[i] = 0.0;
      for (j=0; j<nprimitive(); j++) {
          precon[i] += coef[i][j] * exps[j];
        }
    }

  // precompute contractions over exponentials with exponent weighting
  double precon_g[MAX_NCON];
  if (g_values || h_values) {
      for (i=0; i<ncontraction(); i++) {
          precon_g[i] = 0.0;
          for (j=0; j<nprimitive(); j++) {
              precon_g[i] += exp[j] * coef[i][j] * exps[j];
            }
          precon_g[i] *= 2.0;
        }
    }

  // precompute contractions over exponentials with exponent^2 weighting
  double precon_h[MAX_NCON];
  if (h_values) {
      for (i=0; i<ncontraction(); i++) {
          precon_h[i] = 0.0;
          for (j=0; j<nprimitive(); j++) {
              precon_h[i] += exp[j] * exp[j] * coef[i][j] * exps[j];
            }
          precon_h[i] *= 4.0;
        }
    }

  // compute the shell values
  int i_basis=0;                // Basis function counter
  if (basis_values) {
      for (i=0; i<ncontraction(); i++) {
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
      for (i=0; i<ncontraction(); i++) {
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

  // compute the hessian of the shell values
  if (h_values) {
      int i_hess=0;                // Basis function counter
      for (i=0; i<ncontraction(); i++) {
          // handle s functions with a special case to speed things up
          if (l[i] == 0) {
              double norm_precon_g = precon_g[i];
              double norm_precon_h = precon_h[i];
              // xx
              h_values[i_hess] = norm_precon_h*xs[2] - norm_precon_g;
              i_hess++;
              // yx
              h_values[i_hess] = norm_precon_h*xs[1]*ys[1];
              i_hess++;
              // yy
              h_values[i_hess] = norm_precon_h*ys[2] - norm_precon_g;
              i_hess++;
              // zx
              h_values[i_hess] = norm_precon_h*zs[1]*xs[1];
              i_hess++;
              // zy
              h_values[i_hess] = norm_precon_h*zs[1]*ys[1];
              i_hess++;
              // zz
              h_values[i_hess] = norm_precon_h*zs[2] - norm_precon_g;
              i_hess++;
            }
          else {
              double *cart_h;
              double tmp_cart_h[6*((MAX_AM+1)*(MAX_AM+2))/2];
              if (!puream[i]) {
                  cart_h = &h_values[i_hess];
                }
              else {
                  cart_h = tmp_cart_h;
                }
              CartesianIter *jp = civec[l[i]];
              CartesianIter& j = *jp;
              int i_cart = 0;
              for (j.start(); j; j.next()) {
                  double pre = precon[i];
                  double pre_g = - precon_g[i];
                  double pre_h = precon_h[i];
                  int a = j.a();
                  int b = j.b();
                  int c = j.c();
                  // xx
                  cart_h[i_cart] = pre_h * xs[a+2]*ys[b]*zs[c]
                                 + pre_g * (a+1) * xs[a]*ys[b]*zs[c];
                  if (a>0) {
                      cart_h[i_cart] += pre_g * a*xs[a]*ys[b]*zs[c];
                      if (a>1) cart_h[i_cart] += pre * a*(a-1)
                                               * xs[a-2]*ys[b]*zs[c];
                    }
                  i_cart++;

                  // yx
                  cart_h[i_cart] = pre_h * xs[a+1]*ys[b+1]*zs[c];
                  if (a>0)
                      cart_h[i_cart] += pre_g * a * xs[a-1]*ys[b+1]*zs[c];
                  if (b>0)
                      cart_h[i_cart] += pre_g * b * xs[a+1]*ys[b-1]*zs[c];
                  if (a>0 && b>0)
                      cart_h[i_cart] += pre * a*b * xs[a-1]*ys[b-1]*zs[c];
                  i_cart++;

                  // yy
                  cart_h[i_cart] = pre_h * xs[a]*ys[b+2]*zs[c]
                                 + pre_g * (b+1) * xs[a]*ys[b]*zs[c];
                  if (b>0) {
                      cart_h[i_cart] += pre_g * b*xs[a]*ys[b]*zs[c];
                      if (b>1) cart_h[i_cart] += pre * b*(b-1)
                                               * xs[a]*ys[b-2]*zs[c];
                    }
                  i_cart++;

                  // zx
                  cart_h[i_cart] = pre_h * xs[a+1]*ys[b]*zs[c+1];
                  if (a>0)
                      cart_h[i_cart] += pre_g * a * xs[a-1]*ys[b]*zs[c+1];
                  if (c>0)
                      cart_h[i_cart] += pre_g * c * xs[a+1]*ys[b]*zs[c-1];
                  if (a>0 && c>0)
                      cart_h[i_cart] += pre * a*c * xs[a-1]*ys[b]*zs[c-1];
                  i_cart++;

                  // zy
                  cart_h[i_cart] = pre_h * xs[a]*ys[b+1]*zs[c+1];
                  if (c>0)
                      cart_h[i_cart] += pre_g * c * xs[a]*ys[b+1]*zs[c-1];
                  if (b>0)
                      cart_h[i_cart] += pre_g * b * xs[a]*ys[b-1]*zs[c+1];
                  if (c>0 && b>0)
                      cart_h[i_cart] += pre * c*b * xs[a]*ys[b-1]*zs[c-1];
                  i_cart++;

                  // zz
                  cart_h[i_cart] = pre_h * xs[a]*ys[b]*zs[c+2]
                                 + pre_g * (c+1) * xs[a]*ys[b]*zs[c];
                  if (c>0) {
                      cart_h[i_cart] += pre_g * c*xs[a]*ys[b]*zs[c];
                      if (c>1) cart_h[i_cart] += pre * c*(c-1)
                                               * xs[a]*ys[b]*zs[c-2];
                    }
                  i_cart++;
                }
              if (puream[i]) {
                  SphericalTransformIter *ti = sivec[l[i]];
                  int n = ti->n();
                  memset(&h_values[i_hess], 0, sizeof(double)*n*6);
                  for (ti->start(); ti->ready(); ti->next()) {
                      double coef = ti->coef();
                      int pi = ti->pureindex();
                      int ci = ti->cartindex();
                      for (int xyz2=0; xyz2<6; xyz2++) {
                          h_values[i_hess + pi*6 + xyz2]
                              += coef * cart_h[ci*6 + xyz2];
                        }
                    }
                  i_hess += 6*n;
                }
              else {
                  i_hess += 3*(l[i]+1)*(l[i]+2);
                }
            }
        }
    }

  return i_basis;
}

int
GaussianShell::test_monobound(double &r, double &bound) const
{
  // compute the maximum angular momentum component of the shell
  // add one since derivatives will be needed
  int maxam = max_am() + 1;

  // check limitations
  if (nprimitive() > MAX_NPRIM || ncontraction() > MAX_NCON || maxam >= MAX_AM) {
      std::ostringstream oss;
      oss << "GaussianShell::gaussshval: limit exceeded:\n"
           << indent
           << scprintf(
               "ncon = %d (%d max) nprim = %d (%d max) maxam = %d (%d max)\n",
               ncontraction(),MAX_NCON,nprimitive(),MAX_NPRIM,maxam,MAX_AM-1);
      ExEnv::out0() << indent << oss.str();
      throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
    }

  // loop variables
  int i,j;

  // precompute powers of r
  double rs[MAX_AM+1];
  rs[0] = 1.0;
  if (maxam>0) {
      rs[1] = r;
    }
  for (i=2; i<=maxam; i++) {
      rs[i] = rs[i-1]*r;
    }

  // precompute r*r
  double r2 = r*r;

  // precompute exponentials
  double exps[MAX_NPRIM];
  for (i=0; i<nprimitive(); i++) {
      exps[i]=::exp(-r2*exp[i]);
    }

  // precompute contractions over exponentials
  double precon[MAX_NCON];
  for (i=0; i<ncontraction(); i++) {
      precon[i] = 0.0;
      for (j=0; j<nprimitive(); j++) {
          // using fabs since we want a monotonically decreasing bound
          precon[i] += fabs(coef[i][j]) * exps[j];
        }
    }

  // precompute contractions over exponentials with exponent weighting
  double precon_w[MAX_NCON];
  for (i=0; i<ncontraction(); i++) {
      precon_w[i] = 0.0;
      for (j=0; j<nprimitive(); j++) {
          precon_w[i] += exp[j] * fabs(coef[i][j]) * exps[j];
        }
    }

  double max_bound = 0.0;
  bound = 0.0;
  for (i=0; i<ncontraction(); i++) {
      // using r^l since r^l >= x^a y^b z^c
      double component_bound = rs[l[i]]*precon[i];
      if (l[i] > 0) {
          double d1 = -2.0*rs[l[i]+1]*precon_w[i];
          double d2 = l[i]*rs[l[i]-1]*precon[i];
          if (d1+d2 > 0) {
              // This bound is no good since the contraction is increasing
              // at this position.  Move r out and return to let the driver
              // call again.
              double rold = r;
              r = sqrt(l[i]*precon[i]/(2.0*precon_w[i]));
              if (r<rold+0.01) r = rold+0.01;
              //ExEnv::outn() << "rejected at " << rold << " trying again at "
              //     << r << endl;
              return 1;
            }
        }
      if (component_bound > max_bound) {
          max_bound = component_bound;
      }
    }

  bound = max_bound;
  return 0;
}

double
GaussianShell::monobound(double r) const
{
  // doesn't work at r <= zero
  if (r<=0.001) r = 0.001;
  double b;
  while (test_monobound(r, b));
  return b;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
