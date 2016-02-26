//
// functional.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#include <cmath>
#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/dft/functional.h>

using namespace sc;

#ifndef HAVE_ISNAN
#define isnan(x) ((x)!=(x))
#endif

#define MIN_DENSITY 1.e-14
#define MIN_GAMMA 1.e-24
#define MIN_SQRTGAMMA 1.e-12
#define MAX_ZETA 1.-1.e-12
#define MIN_ZETA -(1.-1.e-12)

///////////////////////////////////////////////////////////////////////////
// utility functions

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

inline static double
norm(double v[3])
{
  return sqrt(dot(v,v));
}

///////////////////////////////////////////////////////////////////////////
// PointInputData

void
PointInputData::compute_derived(int spin_polarized,
                                int need_gradient,
                                int need_hessian)
{
  a.rho_13 = pow(a.rho, 1.0/3.0);
  if (need_gradient) {
      a.gamma = dot(a.del_rho,a.del_rho);
    }
  if (need_hessian) {
      a.lap_rho = a.hes_rho[XX] + a.hes_rho[YY] + a.hes_rho[ZZ];
    }


  if (spin_polarized) {
      b.rho_13 = pow(b.rho, 1.0/3.0);
      if (need_gradient) {
          b.gamma = dot(b.del_rho,b.del_rho);
          gamma_ab = a.del_rho[0]*b.del_rho[0]
                   + a.del_rho[1]*b.del_rho[1] 
                   + a.del_rho[2]*b.del_rho[2];
        }
      if (need_hessian) {
          b.lap_rho = b.hes_rho[XX] + b.hes_rho[YY] + b.hes_rho[ZZ];
        }
    }
  else {
      b = a;
      if (need_gradient) {
          gamma_ab = a.gamma;
        }
    }
}


///////////////////////////////////////////////////////////////////////////
// DenFunctional

static ClassDesc DenFunctional_cd(
  typeid(DenFunctional),"DenFunctional",1,"public SavableState",
  0, 0, 0);

DenFunctional::DenFunctional(StateIn& s):
  SavableState(s)
{
  s.get(a0_);
  s.get(spin_polarized_);
  s.get(compute_potential_);
}

DenFunctional::DenFunctional()
{
  a0_ = 0;
  spin_polarized_ = 0;
  compute_potential_ = 0;
}

DenFunctional::DenFunctional(const Ref<KeyVal>& keyval)
{
  // a0 is usually zero, except for ACM functionals. 
  a0_ = keyval->doublevalue("a0");
  spin_polarized_ = 0;
  compute_potential_ = 0;
}

DenFunctional::~DenFunctional()
{
}

void
DenFunctional::save_data_state(StateOut& s)
{
  s.put(a0_);
  s.put(spin_polarized_);
  s.put(compute_potential_);
}

double
DenFunctional::a0() const
{
  return a0_;
}

int
DenFunctional::need_density_gradient()
{
  return 0;
}

int
DenFunctional::need_density_hessian()
{
  return 0;
}

void
DenFunctional::set_spin_polarized(int i)
{
  spin_polarized_ = i;
}

void
DenFunctional::set_compute_potential(int i)
{
  compute_potential_ = i;
}

void
DenFunctional::gradient(const PointInputData& id, PointOutputData& od,
                        double *grad_f, int acenter,
                        GaussianBasisSet *basis,
                        const double *dmat_a, const double *dmat_b,
                        int ncontrib, const int *contrib,
                        int ncontrib_bf, const int *contrib_bf,
                        const double *bs, const double *bsg,
                        const double *bsh)
{
  int need_gamma_terms = need_density_gradient();
  point(id, od);
  memset(grad_f, 0, sizeof(double)*basis->ncenter()*3);
#if 0
  ExEnv::outn() << scprintf("gradient: rho_a= %12.8f rho_b= %12.8f need_gamma = %d",
                   id.a.rho, id.b.rho, need_gamma_terms) << std::endl;
  ExEnv::outn() << scprintf("  gamma_aa= %12.8f gamma_bb= %12.8f gamma_ab= % 12.8f",
                   id.a.gamma, id.b.gamma, id.gamma_ab) << std::endl;
  ExEnv::outn() << scprintf("  df_drho_a= % 12.8f df_drho_b= % 12.8f",
                   od.df_drho_a, od.df_drho_b) << std::endl;
  ExEnv::outn() << scprintf("  df_dg_aa= % 12.8f df_dg_bb= % 12.8f df_dg_ab= % 12.8f",
                   od.df_dgamma_aa,od.df_dgamma_bb,od.df_dgamma_ab) << std::endl;
#endif

  if (need_gamma_terms) {
      double drhoa = od.df_drho_a;
      double drhob = od.df_drho_b;
      for (int nu=0; nu<ncontrib_bf; nu++) {
          int nut = contrib_bf[nu];
          int nuatom = basis->shell_to_center(basis->function_to_shell(nut));
          double dfa_phi_nu = drhoa * bs[nu];
          double dfb_phi_nu = drhob * bs[nu];
          for (int mu=0; mu<ncontrib_bf; mu++) {
              int mut = contrib_bf[mu];
              int muatom
                  = basis->shell_to_center(basis->function_to_shell(mut));
              if (muatom!=acenter) {
                  int nutmut
                      = (nut>mut?((nut*(nut+1))/2+mut):((mut*(mut+1))/2+nut));
                  double rho_a = dmat_a[nutmut];
                  double rho_b = dmat_b[nutmut];
                  int ixyz;
                  for (ixyz=0; ixyz<3; ixyz++) {
                      double contrib = -2.0*bsg[mu*3+ixyz]
                                     * (rho_a*dfa_phi_nu + rho_b*dfb_phi_nu);
#define hoff(i,j) ((j)<(i)?((i)*((i)+1))/2+(j):((j)*((j)+1))/2+(i))
                      // gamma_aa contrib
                      if (need_gamma_terms) {
                          contrib
                              += 4.0 * od.df_dgamma_aa * rho_a
                              * ( - bsg[mu*3+ixyz]
                                  * ( id.a.del_rho[0]*bsg[nu*3+0]
                                      +id.a.del_rho[1]*bsg[nu*3+1]
                                      +id.a.del_rho[2]*bsg[nu*3+2])
                                  - bs[nu]
                                  * ( id.a.del_rho[0]*bsh[mu*6+hoff(0,ixyz)]
                                      +id.a.del_rho[1]*bsh[mu*6+hoff(1,ixyz)]
                                      +id.a.del_rho[2]*bsh[mu*6+hoff(2,ixyz)]
                                      )
                                  );
                        }
                      // gamma_ab contrib
                      if (need_gamma_terms) {
                          contrib
                              += 2.0 * od.df_dgamma_ab * rho_a
                              * ( - bsg[mu*3+ixyz]
                                  * ( id.b.del_rho[0]*bsg[nu*3+0]
                                      +id.b.del_rho[1]*bsg[nu*3+1]
                                      +id.b.del_rho[2]*bsg[nu*3+2])
                                  - bs[nu]
                                  * ( id.b.del_rho[0]*bsh[mu*6+hoff(0,ixyz)]
                                      +id.b.del_rho[1]*bsh[mu*6+hoff(1,ixyz)]
                                      +id.b.del_rho[2]*bsh[mu*6+hoff(2,ixyz)]
                                      )
                                  );
                          contrib
                              += 2.0 * od.df_dgamma_ab * rho_b
                              * ( - bsg[mu*3+ixyz]
                                  * ( id.a.del_rho[0]*bsg[nu*3+0]
                                      +id.a.del_rho[1]*bsg[nu*3+1]
                                      +id.a.del_rho[2]*bsg[nu*3+2])
                                  - bs[nu]
                                  * ( id.a.del_rho[0]*bsh[mu*6+hoff(0,ixyz)]
                                      +id.a.del_rho[1]*bsh[mu*6+hoff(1,ixyz)]
                                      +id.a.del_rho[2]*bsh[mu*6+hoff(2,ixyz)]
                                      )
                                  );
                        }
                      // gamma_bb contrib
                      if (need_gamma_terms) {
                          contrib
                              += 4.0 * od.df_dgamma_bb * rho_b
                              * ( - bsg[mu*3+ixyz]
                                  * ( id.b.del_rho[0]*bsg[nu*3+0]
                                      +id.b.del_rho[1]*bsg[nu*3+1]
                                      +id.b.del_rho[2]*bsg[nu*3+2])
                                  - bs[nu]
                                  * ( id.b.del_rho[0]*bsh[mu*6+hoff(0,ixyz)]
                                      +id.b.del_rho[1]*bsh[mu*6+hoff(1,ixyz)]
                                      +id.b.del_rho[2]*bsh[mu*6+hoff(2,ixyz)]
                                      )
                                  );
                        }
                      grad_f[3*muatom+ixyz] += contrib;
                      grad_f[3*acenter+ixyz] -= contrib;
                    }
                }
            }
        }
    }
  else {
      double drhoa = od.df_drho_a;
      double drhob = od.df_drho_b;
      for (int nu=0; nu<ncontrib_bf; nu++) {
          int nut = contrib_bf[nu];
          int nuatom = basis->shell_to_center(basis->function_to_shell(nut));
          double dfa_phi_nu = drhoa * bs[nu];
          double dfb_phi_nu = drhob * bs[nu];
          for (int mu=0; mu<ncontrib_bf; mu++) {
              int mut = contrib_bf[mu];
              int muatom
                  = basis->shell_to_center(basis->function_to_shell(mut));
              if (muatom!=acenter) {
                  int nutmut
                      = (nut>mut?((nut*(nut+1))/2+mut):((mut*(mut+1))/2+nut));
                  double rho_a = dmat_a[nutmut];
                  double rho_b = dmat_b[nutmut];
                  int ixyz;
                  for (ixyz=0; ixyz<3; ixyz++) {
//                       std::cout << "bsg[mu*3+ixyz] = " << bsg[mu*3+ixyz]
//                                 << std::endl;
//                       std::cout << "rho_a = " << rho_a
//                                 << std::endl;
//                       std::cout << "rho_b = " << rho_b
//                                 << std::endl;
//                       std::cout << "dfa_phi_nu = " << dfa_phi_nu
//                                 << std::endl;
//                       std::cout << "dfb_phi_nu = " << dfb_phi_nu
//                                 << std::endl;
                      double contrib = -2.0*bsg[mu*3+ixyz]
                                     * (rho_a*dfa_phi_nu + rho_b*dfb_phi_nu);
                      grad_f[3*muatom+ixyz] += contrib;
                      grad_f[3*acenter+ixyz] -= contrib;
                    }
                }
            }
        }
    }
}

static void
compute_derived_without_del_rho(PointInputData &id,
                                int spin_polarized,
                                int need_gradient,
                                int need_hessian)
{
  id.a.rho_13 = pow(id.a.rho, 1.0/3.0);

  if (!spin_polarized) {
      id.b.rho = id.a.rho;
      id.b.gamma = id.a.gamma;
      id.gamma_ab = id.a.gamma;
    }

  id.b.rho_13 = pow(id.b.rho, 1.0/3.0);

  // These quantities are not initialized and so are set to nonsense
  // values.  Since no functional of which I am uses del_rho directly, this
  // should not affect the results.
  for (int i=0; i<3; i++) {
      id.a.del_rho[i] = 1e99;
      id.b.del_rho[i] = 1e99;
    }
  for (int i=0; i<6; i++) {
      id.a.hes_rho[i] = 1e99;
      id.b.hes_rho[i] = 1e99;
    }
}

void
DenFunctional::do_fd_point(PointInputData&id,
                           double&in,double&out,
                           double lower_bound, double upper_bound)
{
  double delta = 0.0000000001;
  PointOutputData tod;
  double insave = in;

  point(id,tod);
  double outsave = tod.energy;

  if (insave-delta>=lower_bound && insave+delta<=upper_bound) {
      in = insave+delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double plus = tod.energy;

      in = insave-delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double minu = tod.energy;
      out = 0.5*(plus-minu)/delta;
    }
  else if (insave+2*delta<=upper_bound) {
      in = insave+delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double plus = tod.energy;

      in = insave+2*delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double plus2 = tod.energy;
      out = 0.5*(4.0*plus-plus2-3.0*outsave)/delta;
    }
  else if (insave-2*delta>=lower_bound) {
      in = insave-delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double minu = tod.energy;

      in = insave-2*delta;
      compute_derived_without_del_rho(id, spin_polarized_,
                                      need_density_gradient(), false);
      point(id,tod);
      double minu2 = tod.energy;
      out = -0.5*(4.0*minu-minu2-3.0*outsave)/delta;
    }
  else {
      // the derivative is not well defined for this case
      out = -135711.;
    }
  in = insave;
  compute_derived_without_del_rho(id, spin_polarized_,
                                  need_density_gradient(), false);
}

void
DenFunctional::fd_point(const PointInputData&id, PointOutputData&od)
{
  PointInputData tid(id);

  // fill in the energy at the initial density values
  point(id,od);

  ExEnv::out0() << scprintf("ra=%7.5f rb=%7.5f gaa=%9.7f gbb=%9.7f gab= % 9.7f",
                   id.a.rho, id.b.rho, id.a.gamma, id.b.gamma, id.gamma_ab)
       << std::endl;

  if (spin_polarized_) {
      double ga = tid.a.gamma;
      double gb = tid.b.gamma;
      double gab = tid.gamma_ab;
      double sga = sqrt(ga);
      double sgb = sqrt(gb);

      double g_a_lbound = -2*gab - gb;
      if (gb > 0 && gab*gab/gb > g_a_lbound) g_a_lbound = gab*gab/gb;
      if (g_a_lbound < 0) g_a_lbound = 0.0;

      double g_b_lbound = -2*gab - ga;
      if (ga > 0 && gab*gab/ga > g_b_lbound) g_b_lbound = gab*gab/ga;
      if (g_b_lbound < 0) g_b_lbound = 0.0;

      double g_ab_lbound = -0.5*(ga+gb);
      if (-sga*sgb > g_ab_lbound) g_ab_lbound = -sga*sgb;
      // if (-sga*sgb < g_ab_lbound) g_ab_lbound = -sga*sgb;
      double g_ab_ubound = sga*sgb;

      do_fd_point(tid, tid.a.rho, od.df_drho_a, 0.0, 10.0);
      do_fd_point(tid, tid.b.rho, od.df_drho_b, 0.0, 10.0);
      do_fd_point(tid, tid.a.gamma, od.df_dgamma_aa, g_a_lbound, 10.0);
      do_fd_point(tid, tid.b.gamma, od.df_dgamma_bb, g_b_lbound, 10.0);
      do_fd_point(tid, tid.gamma_ab, od.df_dgamma_ab, g_ab_lbound, g_ab_ubound);
    }
  else {
      double delta = 0.0000000001;
      PointOutputData tod;
      double plus, plus2;

      // Do the rho displacements.
      // Use two positive displacements to obtain good accuracy
      // and to avoid negative values.
      // Since the alpha displacement is copied to beta, the
      // delta must be divided by 2 to get the desired displacement
      // in rho = rho_alpha + rho_beta.
      tid.a.rho += (delta/2.0);
      compute_derived_without_del_rho(tid, spin_polarized_,
                                      need_density_gradient(), false);
      point(tid,tod);
      plus = tod.energy;
      tid.a.rho = id.a.rho;
      tid.a.rho += 2.0*(delta/2.0);
      compute_derived_without_del_rho(tid, spin_polarized_,
                                      need_density_gradient(), false);
      point(tid,tod);
      plus2 = tod.energy;
      double df_drho = 0.5*(4.0*plus-plus2-3.0*od.energy)/delta;
      od.df_drho_a = df_drho;
      od.df_drho_b = df_drho;
      tid.a.rho = id.a.rho;

      // Do the gamma displacements.
      // Use two positive displacements to obtain good accuracy
      // and to avoid negative values.
      // Since the alpha, alpha displacement is copied to beta, beta and
      // alpha, beta, the delta must be divided by 4 to get the desired
      // displacement in gamma = gamma_alpha,alpha + gamma_beta,beta +
      // 2 gamma_alpha,beta.
      tid.a.gamma += (delta/4.0);
      compute_derived_without_del_rho(tid, spin_polarized_,
                                      need_density_gradient(), false);
      point(tid,tod);
      plus = tod.energy;
      tid.a.gamma = id.a.gamma;
      tid.a.gamma += 2.0*(delta/4.0);
      compute_derived_without_del_rho(tid, spin_polarized_,
                                      need_density_gradient(), false);
      point(tid,tod);
      plus2 = tod.energy;
      tid.a.gamma = id.a.gamma;
      double df_dgamma = 0.5*(4.0*plus-plus2-3.0*od.energy)/delta;
      od.df_dgamma_aa =    df_dgamma;
      od.df_dgamma_bb =    df_dgamma;
      od.df_dgamma_ab = 2.*df_dgamma;
      compute_derived_without_del_rho(tid, spin_polarized_,
                                      need_density_gradient(), false);
    }
}

static int
check(const char *name, double fd, double an, const char *class_name)
{
  double threshold = 1.e-4;
  double tolerance = 5.e-2;
  double err = fabs(fd - an);
  // -135711. flags an undefined FD
  if (fd == -135711.) {
      ExEnv::out0() << scprintf("%20s: fd = %12s an = % 12.8f",
                                name, "undefined", an)
                    << std::endl;
    }
  else if ((fabs(an) > threshold && err/fabs(an) > tolerance)
      || ((fabs(an) <= threshold) && err > tolerance)
#ifdef HAVE_ISNAN
      || isnan(an)
#endif
      ) {
      ExEnv::out0() << scprintf("\033[31mERROR: %13s: fd = % 12.8f an = % 12.8f (%s)\033[30m",
                       name, fd, an, class_name)
           << std::endl;
      return 1;
    }
  else {
      ExEnv::out0() << scprintf("%20s: fd = % 12.8f an = % 12.8f", name, fd, an)
                    << std::endl;
    }
  return 0;
}

int
DenFunctional::test(const PointInputData &id)
{
  PointOutputData fd_od;
  fd_point(id,fd_od);
  PointOutputData an_od;
  point(id,an_od);
  int r = 0;
  ExEnv::out0() << scprintf("%20s  f = % 12.8f", "", an_od.energy) << std::endl;
  r+=check("df_drho_a", fd_od.df_drho_a, an_od.df_drho_a, class_name());
  r+=check("df_drho_b", fd_od.df_drho_b, an_od.df_drho_b, class_name());
  r+=check("df_dgamma_aa",fd_od.df_dgamma_aa,an_od.df_dgamma_aa,class_name());
  r+=check("df_dgamma_ab",fd_od.df_dgamma_ab,an_od.df_dgamma_ab,class_name());
  r+=check("df_dgamma_bb",fd_od.df_dgamma_bb,an_od.df_dgamma_bb,class_name());
  return r;
}

int
DenFunctional::test()
{
  int i, j, k, l, m;
  set_compute_potential(1);
  SCVector3 r = 0.0;
  PointInputData id(r);

  for (i=0; i<6; i++) {
      id.a.hes_rho[i] = 0.0;
      id.b.hes_rho[i] = 0.0;
    }
  id.a.lap_rho = 0.0;
  id.b.lap_rho = 0.0;

  std::vector<double> testrho(6);
  testrho[0] = 0.000;
  testrho[1] = 0.001;
  testrho[2] = 0.010;
  testrho[3] = 0.100;
  testrho[4] = 0.500;
  testrho[5] = 1.000;
  std::vector<std::vector<double> > testdelrho(3);
  for (int i=0; i<testdelrho.size(); i++)
      testdelrho[i].resize(3);
  testdelrho[0][0] = 0.000;
  testdelrho[0][1] = 0.000;
  testdelrho[0][2] = 0.000;
  testdelrho[1][0] = 0.001;
  testdelrho[1][1] = 0.000;
  testdelrho[1][2] = 0.000;
  testdelrho[2][0] = 0.000;
  testdelrho[2][1] = 0.500;
  testdelrho[2][2] = 0.000;

  int ret = 0;

  ExEnv::out0() << "Testing with rho_a == rho_b" << std::endl;
  set_spin_polarized(0);
  for (i=0; i < testrho.size(); i++) {
      if (testrho[i] == 0.0) continue;
      id.a.rho=testrho[i];
      id.b.rho=testrho[i];
      for (j=0; j < testdelrho.size(); j++) {
          for (int k=0; k<3; k++) id.a.del_rho[k] = testdelrho[j][k];
          id.compute_derived(0, need_density_gradient(), false);
          // constrain gamma to be more or less physically reasonable
          if (sqrt(id.a.gamma) > 1.e3*id.a.rho) continue;
          ExEnv::out0() << "testing rho = " << id.a.rho
                        << " delrho = {" << id.a.del_rho[0]
                        << ", " << id.a.del_rho[1]
                        << ", " << id.a.del_rho[2]
                        << "} gamma = " << id.a.gamma
                        << std::endl;
          ret += test(id);
        }
    }

  try {
      set_spin_polarized(1);
      ExEnv::out0() << "Testing with rho_a != rho_b" << std::endl;
      for (i=0; i<testrho.size(); i++) {
          id.a.rho=testrho[i];
          for (j=0; j<testrho.size(); j++) {
              id.b.rho=testrho[j];
              if (testrho[i]+testrho[j] == 0.0) continue;
              for (k=0; k<testdelrho.size(); k++) {
                  for (int xyz=0;xyz<3;xyz++)
                      id.a.del_rho[xyz] = testdelrho[k][xyz];
                  for (l=0; l < testdelrho.size(); l++) {
                      for (int xyz=0;xyz<3;xyz++)
                          id.b.del_rho[xyz] = testdelrho[l][xyz];
                      id.compute_derived(1, need_density_gradient(), false);
                      // constrain gamma to be more or less physically reasonable
                      if (sqrt(id.a.gamma) > 2.*id.a.rho) continue;
                      if (sqrt(id.b.gamma) > 2.*id.b.rho) continue;
                      ret += test(id);
                    }
                }
            }
        }
    }
  catch(std::exception&e) {
      ExEnv::out0() << indent
                    << "Caught an exception when testing with rho_a != rho_b"
                    << std::endl
                    << indent
                    << e.what()
                    << std::endl
                    << indent
                    << "skipping test"
                    << std::endl;

    }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////
// NElFunctional

static ClassDesc NElFunctional_cd(
  typeid(NElFunctional),"NElFunctional",1,"public DenFunctional",
  0, create<NElFunctional>, create<NElFunctional>);

NElFunctional::NElFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

NElFunctional::NElFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
}

NElFunctional::~NElFunctional()
{
}

void
NElFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
}

void
NElFunctional::point(const PointInputData &id,
                     PointOutputData &od)
{
  od.zero();
  od.energy = id.a.rho + id.b.rho;
}

/////////////////////////////////////////////////////////////////////////////
// SumDenFunctional

static ClassDesc SumDenFunctional_cd(
  typeid(SumDenFunctional),"SumDenFunctional",1,"public DenFunctional",
  0, create<SumDenFunctional>, create<SumDenFunctional>);

SumDenFunctional::SumDenFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s),
  n_(0),
  funcs_(0),
  coefs_(0)
{
  s.get(n_);
  if (n_) {
      s.get(coefs_);
      funcs_ = new Ref<DenFunctional>[n_];
      for (int i=0; i < n_; i++)
          funcs_[i] << SavableState::restore_state(s);
    }
}

SumDenFunctional::SumDenFunctional() :
  n_(0),
  funcs_(0),
  coefs_(0)
{
}

SumDenFunctional::SumDenFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval),
  n_(0),
  funcs_(0),
  coefs_(0)
{
  int ncoef = keyval->count("coefs");
  int nfunc = keyval->count("funcs");
  if (ncoef != nfunc && ncoef != 0) {
      ExEnv::out0() << "SumDenFunctional: number of coefs and funcs differ" << std::endl;
      abort();
    }
  
  n_ = nfunc;
  coefs_ = new double[n_];
  funcs_ = new Ref<DenFunctional>[n_];
  for (int i=0; i < n_; i++) {
      if (ncoef)
          coefs_[i] = keyval->doublevalue("coefs", i);
      else
          coefs_[i] = 1.0;
      funcs_[i] << keyval->describedclassvalue("funcs", i);
    }
}

SumDenFunctional::~SumDenFunctional()
{
  if (n_) {
      for (int i=0; i < n_; i++) funcs_[i] = 0; // just in case
      delete[] funcs_;
      delete[] coefs_;
    }
  n_=0;
  funcs_=0;
  coefs_=0;
}

void
SumDenFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(n_);
  if (n_) {
      s.put(coefs_, n_);
      for (int i=0; i < n_; i++) 
          SavableState::save_state(funcs_[i].pointer(),s);
    }
}

double
SumDenFunctional::a0() const
{
  double eff_a0 = a0_;
  for (int i=0; i < n_; i++) {
      eff_a0 += coefs_[i] * funcs_[i]->a0();
    }
  return eff_a0;
}

int
SumDenFunctional::need_density_gradient()
{
  for (int i=0; i < n_; i++)
      if (funcs_[i]->need_density_gradient())
          return 1;

  return 0;
}

void
SumDenFunctional::set_spin_polarized(int p)
{
  spin_polarized_ = p;
  for (int i=0; i < n_; i++)
      funcs_[i]->set_spin_polarized(p);
}

void
SumDenFunctional::set_compute_potential(int val)
{
  compute_potential_ = val;
  for (int i=0; i < n_; i++)
      funcs_[i]->set_compute_potential(val);
}

void
SumDenFunctional::point(const PointInputData &id,
                        PointOutputData &od)
{
  od.zero();
  PointOutputData tmpod;
  for (int i=0; i < n_; i++) {
      funcs_[i]->point(id, tmpod);
      
      od.energy += coefs_[i] * tmpod.energy;
      if (compute_potential_) {
          od.df_drho_a += coefs_[i] * tmpod.df_drho_a;
          od.df_drho_b += coefs_[i] * tmpod.df_drho_b;
          od.df_dgamma_aa += coefs_[i] * tmpod.df_dgamma_aa;
          od.df_dgamma_ab += coefs_[i] * tmpod.df_dgamma_ab;
          od.df_dgamma_bb += coefs_[i] * tmpod.df_dgamma_bb;
        }
    }
}

void
SumDenFunctional::print(std::ostream& o) const
{
  o
    << indent << "Sum of Functionals:" << std::endl;
  o << incindent;
  o << indent << scprintf("%+18.16f Hartree-Fock Exchange",a0_) << std::endl;
  for (int i=0; i<n_; i++) {
      o << indent << scprintf("%+18.16f",coefs_[i]) << std::endl;
      o << incindent;
      funcs_[i]->print(o);
      o << decindent;
    }
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////
// StdDenFunctional

static ClassDesc StdDenFunctional_cd(
  typeid(StdDenFunctional),"StdDenFunctional",1,"public SumDenFunctional",
  0, create<StdDenFunctional>, create<StdDenFunctional>);

StdDenFunctional::StdDenFunctional(StateIn& s):
  SavableState(s),
  SumDenFunctional(s)
{
  s.get(name_);
}

StdDenFunctional::StdDenFunctional()
{
}

void
StdDenFunctional::init_arrays(int n)
{
  n_ = n;
  funcs_ = new Ref<DenFunctional>[n_];
  coefs_ = new double[n_];
  for (int i=0; i<n_; i++) coefs_[i] = 1.0;
}

StdDenFunctional::StdDenFunctional(const Ref<KeyVal>& keyval)
{
  name_ = keyval->stringvalue("name");
  if (!name_.empty()) {
      if (name_ == "HFK") {
          n_ = 0;
          a0_ = 1.0;
        }
      else if (name_ == "XALPHA") {
          init_arrays(1);
          funcs_[0] = new XalphaFunctional;
        }
      else if (name_ == "HFS") {
          init_arrays(1);
          funcs_[0] = new SlaterXFunctional;
        }
      else if (name_ == "HFB") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
        }
      else if (name_ == "HFG96") {
          init_arrays(1);
          funcs_[0] = new G96XFunctional;
        }
      else if (name_ == "G96LYP") {
          init_arrays(2);
          funcs_[0] = new G96XFunctional;
          funcs_[1] = new LYPCFunctional;
        }
      else if (name_ == "BLYP") {
          init_arrays(3);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new LYPCFunctional;
        }
      else if (name_ == "SVWN1") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN1LCFunctional;
        }
      else if (name_ == "SVWN1RPA") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN1LCFunctional(1);
        }
      else if (name_ == "SVWN2") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN2LCFunctional;
        }
      else if (name_ == "SVWN3") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN3LCFunctional;
        }
      else if (name_ == "SVWN4") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN4LCFunctional;
        }
      else if (name_ == "SVWN5") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN5LCFunctional;
        }
      else if (name_ == "SPZ81") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new PZ81LCFunctional;
        }
      else if (name_ == "SPW92") {
          init_arrays(2);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new PW92LCFunctional;
        }
      else if (name_ == "BPW91") {
          init_arrays(3);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new PW91CFunctional;
        }
      else if (name_ == "BP86") {
          init_arrays(4);
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new P86CFunctional;
          funcs_[3] = new PZ81LCFunctional;
        }
      else if (name_ == "B3LYP") {
          init_arrays(4);
          a0_ = 0.2;
          coefs_[0] = 0.8;
          coefs_[1] = 0.72;
          coefs_[2] = 0.19;
          coefs_[3] = 0.81;
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new VWN1LCFunctional(1);
          funcs_[3] = new LYPCFunctional;
        }
      else if (name_ == "KMLYP") {
          init_arrays(3);
          a0_ = 0.557;
          coefs_[0] = 0.443;
          coefs_[1] = 0.552;
          coefs_[2] = 0.448;
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new VWN1LCFunctional(1);
          funcs_[2] = new LYPCFunctional;
        }
      else if (name_ == "B3PW91") {
          init_arrays(4);
          a0_ = 0.2;
          coefs_[0] = 0.8;
          coefs_[1] = 0.72;
          coefs_[2] = 0.81;
          coefs_[3] = 0.19;
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new PW91CFunctional;
          funcs_[3] = new PW92LCFunctional;
        }
      else if (name_ == "B3P86") {
          init_arrays(4);
          a0_ = 0.2;
          coefs_[0] = 0.8;
          coefs_[1] = 0.72;
          coefs_[2] = 0.81;
          coefs_[3] = 1.0;
          funcs_[0] = new SlaterXFunctional;
          funcs_[1] = new Becke88XFunctional;
          funcs_[2] = new P86CFunctional;
          funcs_[3] = new VWN1LCFunctional(1);
        }
      else if (name_ == "PBE") {
          init_arrays(2);
          funcs_[0] = new PBEXFunctional;
          funcs_[1] = new PBECFunctional;
        }
      else if (name_ == "PW91") {
          init_arrays(2);
          funcs_[0] = new PW91XFunctional;
          funcs_[1] = new PW91CFunctional;
        }
      else if (name_ == "mPW(PW91)PW91") {
          init_arrays(2);
          funcs_[0] = new mPW91XFunctional(mPW91XFunctional::PW91);
          funcs_[1] = new PW91CFunctional;
        }
      else if (name_ == "mPWPW91") {
          init_arrays(2);
          funcs_[0] = new mPW91XFunctional(mPW91XFunctional::mPW91);
          funcs_[1] = new PW91CFunctional;
        }
      else if (name_ == "mPW1PW91") {
          init_arrays(2);
          a0_ = 0.16;
          coefs_[0] = 0.84;
          coefs_[1] = 1.0;
          funcs_[0] = new mPW91XFunctional(mPW91XFunctional::mPW91);
          funcs_[1] = new PW91CFunctional;
        }
      else {
          ExEnv::out0() << "StdDenFunctional: bad name: " << name_ << std::endl;
          abort();
        }
    }
}

StdDenFunctional::~StdDenFunctional()
{
}

void
StdDenFunctional::save_data_state(StateOut& s)
{
  SumDenFunctional::save_data_state(s);
  s.put(name_);
}

void
StdDenFunctional::print(std::ostream& o) const
{
  o
    << indent << "Standard Density Functional: " << name_ << std::endl;
  SumDenFunctional::print(o);
}

/////////////////////////////////////////////////////////////////////////////
// LSDACFunctional: All local correlation functionals inherit from this class.
// Coded by Matt Leininger
static ClassDesc LSDACFunctional_cd(
  typeid(LSDACFunctional),"LSDACFunctional",1,"public DenFunctional",
  0, 0, 0);

LSDACFunctional::LSDACFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

LSDACFunctional::LSDACFunctional()
{
}

LSDACFunctional::LSDACFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
}

LSDACFunctional::~LSDACFunctional()
{
}

void
LSDACFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
}

void
LSDACFunctional::point(const PointInputData &id,
                       PointOutputData &od)
{
  double junk_1, junk_2, junk_3;
  point_lc(id, od, junk_1, junk_2, junk_3);
}

/////////////////////////////////////////////////////////////////////////////
// SlaterXFunctional

static ClassDesc SlaterXFunctional_cd(
  typeid(SlaterXFunctional),"SlaterXFunctional",1,"public DenFunctional",
  0, create<SlaterXFunctional>, create<SlaterXFunctional>);

SlaterXFunctional::SlaterXFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

SlaterXFunctional::SlaterXFunctional()
{
}

SlaterXFunctional::SlaterXFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
}

SlaterXFunctional::~SlaterXFunctional()
{
}

void
SlaterXFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
}

void
SlaterXFunctional::point(const PointInputData &id,
                       PointOutputData &od)
{
  const double mcx2rthird = -0.9305257363491; // -1.5*(3/4pi)^1/3
  const double dmcx2rthird = -1.2407009817988; // 2*(3/4pi)^1/3
  od.zero();

  if (!spin_polarized_) {
      od.energy = mcx2rthird * 2.0 * id.a.rho * id.a.rho_13;
      if (compute_potential_) {
          od.df_drho_a = dmcx2rthird * id.a.rho_13;
          od.df_drho_b = od.df_drho_a;
          }
    }
  else {
      od.energy = mcx2rthird
                * (id.a.rho * id.a.rho_13 + id.b.rho * id.b.rho_13);
      if (compute_potential_) {
          od.df_drho_a = dmcx2rthird * id.a.rho_13;
          od.df_drho_b = dmcx2rthird * id.b.rho_13;
          }
    }
}

/////////////////////////////////////////////////////////////////////////////
// PW92LCFunctional
// Coded by Matt Leininger
static ClassDesc PW92LCFunctional_cd(
  typeid(PW92LCFunctional),"PW92LCFunctional",1,"public LSDACFunctional",
  0, create<PW92LCFunctional>, create<PW92LCFunctional>);

PW92LCFunctional::PW92LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

PW92LCFunctional::PW92LCFunctional()
{
}

PW92LCFunctional::PW92LCFunctional(const Ref<KeyVal>& keyval):
  LSDACFunctional(keyval)
{
}

PW92LCFunctional::~PW92LCFunctional()
{
}

void
PW92LCFunctional::save_data_state(StateOut& s)
{
  LSDACFunctional::save_data_state(s);
}

double
PW92LCFunctional::F(double x, double A, double alpha_1, double beta_1, double beta_2,
                    double beta_3, double beta_4, double p)
{
  double x2 = x*x; // r_s
  double denom = 2.*A*( beta_1 * x + beta_2 * x2 + beta_3 * x2*x + beta_4 * pow(x2,p+1.));
  double res = -2.*A*(1. + alpha_1*x2)*log(1.+ 1./denom);

  return res;
}

double
PW92LCFunctional::dFdr_s(double x, double A, double alpha_1, double beta_1, double beta_2,
                    double beta_3, double beta_4, double p)
{
  double x2 = x*x; // r_s
  double Q_0 = -2.*A*(1. + alpha_1*x2);
  double Q_1 =  2.*A*(beta_1 * x + beta_2 * x2 + beta_3*x*x2 + beta_4 * pow(x2,p+1.));
  double Q_1prime = A * 
           ( beta_1 * 1./x + 2.*beta_2 + 3.*beta_3*x + 2.*(p+1.)*beta_4*pow(x2,p));
  double res = -2.*A*alpha_1*log(1. + 1./Q_1) - Q_0*Q_1prime/(Q_1*Q_1 + Q_1);

  return res;
}

void
PW92LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                           double &ec_local, double &decrs, double &deczeta)
{
  od.zero();
  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  const double sixth = 1./6.;
  const double four_thirds = 4./3.;
  const double one_third = 1./3.;
  const double two_thirds = 2./3.;

  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;
  
  double epc    = F(x, 0.0310907,  0.21370, 7.5957,  3.5876, 1.6382,  0.49294, 1.00);
  double efc    = F(x, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662,  0.62517, 1.00);
  double alphac = F(x, 0.0168869, 0.11125, 10.357,  3.6231, 0.88026, 0.49671, 1.00);
     
  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_ec = -alphac * f / fpp0 * (1. - zeta4) + (efc - epc) * f * zeta4;
  double ec = epc + delta_ec;

  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
    if (!spin_polarized_) {
      double depc_dr_s0 = 
             dFdr_s(x, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, 1.00);
      double dec_dr_s = depc_dr_s0;
      od.df_drho_a = od.df_drho_b = ec - (rs/3.)*dec_dr_s;
      decrs = dec_dr_s;
      deczeta = 0.;
    }
    else {
      double zeta3 = zeta2*zeta;
      double depc_dr_s0 = dFdr_s(x, 0.0310907, 0.21370, 7.5957,  
                                    3.5876, 1.6382,  0.49294, 1.00);
      double defc_dr_s1 = dFdr_s(x, 0.01554535, 0.20548, 14.1189, 
                                    6.1977, 3.3662,  0.62517, 1.00);
      double dalphac_dr_s = dFdr_s(x, 0.0168869, 0.11125, 10.357,  
                                    3.6231, 0.88026, 0.49671, 1.00);
      double dec_dr_s = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
                        + -dalphac_dr_s * f / fpp0 * (1 - zeta4);
      double fp = two_thirds * (pow((1+zeta),one_third) 
            - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      double dec_dzeta = 4.* zeta3 * f * (efc - epc - (-alphac/fpp0))
              + fp * (zeta4 * (efc - epc) + (1-zeta4)*(-alphac/fpp0));
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
      } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// PZ81LCFunctional
// Coded by Matt Leininger
// Used in P86 correlation functional
// J. P. Perdew and A. Zunger, Phys. Rev. B, 23, 5048, 1981.
// C. W. Murray, N. C. Handy, G. J. Laming, Mol. Phys., 78, 997, 1993.
//
static ClassDesc PZ81LCFunctional_cd(
  typeid(PZ81LCFunctional),"PZ81LCFunctional",1,"public LSDACFunctional",
  0, create<PZ81LCFunctional>, create<PZ81LCFunctional>);

PZ81LCFunctional::PZ81LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

PZ81LCFunctional::PZ81LCFunctional()
{
}

PZ81LCFunctional::PZ81LCFunctional(const Ref<KeyVal>& keyval):
  LSDACFunctional(keyval)
{
}

PZ81LCFunctional::~PZ81LCFunctional()
{
}

void
PZ81LCFunctional::save_data_state(StateOut& s)
{
  LSDACFunctional::save_data_state(s);
}

double
PZ81LCFunctional::Fec_rsgt1(double rs, double beta1, double beta2, double gamma)
{
  double sqrt_rs = sqrt(rs);
  double res = gamma / (1. + beta1*sqrt_rs + beta2*rs);
      
  return res;
}

double
PZ81LCFunctional::dFec_rsgt1_drho(double rs, double beta1, double beta2, double gamma,
                                  double &dec_drs)
{
  double ec = Fec_rsgt1(rs, beta1, beta2, gamma);
  double sqrt_rs = sqrt(rs);
  // double numer = 1.+ 7./6.*beta1*sqrt_rs + 4./3.*beta2*rs;
  double denom = 1. + beta1*sqrt_rs + beta2*rs;
  dec_drs = -ec/denom * (beta1/(2.*sqrt_rs) + beta2);
  // double res = ec * numer / denom;
  double res = (ec - rs/3.*dec_drs);
  
  return res;
}

double
PZ81LCFunctional::Fec_rslt1(double rs, double A, double B, double C, double D)
{
  double lnrs = log(rs);
  double res = A*lnrs + B + C*rs*lnrs + D*rs;
      
  return res;
}

double
PZ81LCFunctional::dFec_rslt1_drho(double rs, double A, double B, double C, double D,
                                  double &dec_drs)
{
  double lnrs = log(rs);
  double res = A*lnrs + B - A/3. + 2./3.*C*rs*lnrs + 1./3.*(2.*D - C)*rs;
  dec_drs = A/rs + C*lnrs + C + D;
  return res;
}


void
PZ81LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();

  const double Au = 0.0311;
  const double Ap = 0.01555;
  const double Bu = -0.048;
  const double Bp = -0.0269;
  const double Cu = 0.0020;
  const double Cp = 0.0007;
  const double Du = -0.0116;
  const double Dp = -0.0048;
  const double beta1u = 1.0529;
  const double beta1p = 1.3981;
  const double beta2u = 0.3334;
  const double beta2p = 0.2611;
  const double gammau = -0.1423;
  const double gammap = -0.0843;
  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double rs = pow(3./(4.*M_PI*rho), (1./3.) );
  double fzeta = ( pow((1.+zeta), (4./3.)) + pow((1.-zeta), (4./3.)) - 2.)
                 / ( pow(2., (4./3.)) - 2. );

  double euc, epc;
  if (rs >= 1.) {
    // Ceperley U
    // euc = Fec_rsgt1(rs, 1.1581, 0.3446, -0.1471);
    // Ceperley P
    // epc = Fec_rsgt1(rs, 1.2520, 0.2567, -0.0790);
    // Ceperley-Adler U
    euc = Fec_rsgt1(rs, beta1u, beta2u, gammau);
    // Ceperley-Adler P
    epc = Fec_rsgt1(rs, beta1p, beta2p, gammap);
    }
  else { // rs < 1.
      // Ceperley U with A_u and B_u
      // euc = Fec_rslt1(rs, 0.0311, -0.048, 0.0014, -0.0108);
      // Ceperley P with A_p and B_p
      // epc = Fec_rslt1(rs, 0.01555, -0.0269, 0.0001, -0.0046);
      // Ceperley-Adler U with A_u and B_u
      euc = Fec_rslt1(rs, Au, Bu, Cu, Du);
      // Ceperley-Adler P with A_p and B_p
      epc = Fec_rslt1(rs, Ap, Bp, Cp, Dp);
    }
  double ec = euc + fzeta*(epc-euc);
  ec_local = ec;
  od.energy = ec * rho;

  if (compute_potential_) {
      double deuc_drs = 0.;
      double depc_drs = 0.;
      double deuc_drho, depc_drho;
      if (rs > 1.) {
          // Ceperley U
          // deuc_drho = dFec_rsgt1_drho(rs, 1.1581, 0.3446, -0.1471, deuc_drs);
          // Ceperley P
          // depc_drho = dFec_rsgt1_drho(rs, 1.2520, 0.2567, -0.0790, depc_drs);
          // Ceperley-Adler U
          deuc_drho = dFec_rsgt1_drho(rs, beta1u, beta2u, gammau, deuc_drs);
          // Ceperley-Adler P
          depc_drho = dFec_rsgt1_drho(rs, beta1p, beta2p, gammap, depc_drs);
         }
      else { // rs < 1.
         // Ceperley U with A_u and B_u
         // deuc_drho = dFec_rslt1_drho(rs, 0.0311, -0.048, 0.0014, -0.0108, deuc_drs);
         // Ceperley P with A_p and B_p
         // depc_drho = dFec_rslt1_drho(rs, 0.01555, -0.0269, 0.0001, -0.0046, depc_drs);
         // Ceperley-Adler U with A_u and B_u
         deuc_drho = dFec_rslt1_drho(rs, Au, Bu, Cu, Du, deuc_drs);
         // Ceperley-Adler P with A_p and B_p
         depc_drho = dFec_rslt1_drho(rs, Ap, Bp, Cp, Dp, depc_drs);
         }
      double dfzeta_dzeta = 4./3.*( pow((1.+zeta), (1./3.)) - pow((1.-zeta), (1./3.)) )
                      / (pow(2., (4./3.)) - 2.);
      decrs = deuc_drs + fzeta*(depc_drs - deuc_drs);
      deczeta = dfzeta_dzeta*(epc - euc);
      od.df_drho_a = deuc_drho + fzeta*(depc_drho - deuc_drho) +
                     (epc - euc)*(1.-zeta)*dfzeta_dzeta;
      od.df_drho_b = deuc_drho + fzeta*(depc_drho - deuc_drho) +
                     (epc - euc)*(-1.-zeta)*dfzeta_dzeta;
           
      } 
}

/////////////////////////////////////////////////////////////////////////////
// VWNLCFunctional
// Coded by Matt Leininger

static ClassDesc VWNLCFunctional_cd(
  typeid(VWNLCFunctional),"VWNLCFunctional",1,"public LSDACFunctional",
  0, create<VWNLCFunctional>, create<VWNLCFunctional>);

VWNLCFunctional::VWNLCFunctional(StateIn& s):
  SavableState(s), LSDACFunctional(s)
{
  s.get(Ap_);
  s.get(Af_);
  s.get(A_alpha_);

  s.get(x0p_mc_);
  s.get(bp_mc_);
  s.get(cp_mc_);
  s.get(x0f_mc_);
  s.get(bf_mc_);
  s.get(cf_mc_);
  s.get(x0_alpha_mc_);
  s.get(b_alpha_mc_);
  s.get(c_alpha_mc_);
  
  s.get(x0p_rpa_);
  s.get(bp_rpa_);
  s.get(cp_rpa_);
  s.get(x0f_rpa_);
  s.get(bf_rpa_);
  s.get(cf_rpa_);
  s.get(x0_alpha_rpa_);
  s.get(b_alpha_rpa_);
  s.get(c_alpha_rpa_);
}

VWNLCFunctional::VWNLCFunctional()
{
  init_constants();
}

VWNLCFunctional::~VWNLCFunctional()
{
}

VWNLCFunctional::VWNLCFunctional(const Ref<KeyVal>& keyval):
  LSDACFunctional(keyval)
{
  init_constants();
  Ap_  = keyval->doublevalue("Ap", KeyValValuedouble(Ap_));
  Af_  = keyval->doublevalue("Af", KeyValValuedouble(Af_));
  A_alpha_  = keyval->doublevalue("A_alpha", KeyValValuedouble(A_alpha_));
 
  x0p_mc_  = keyval->doublevalue("x0p_mc", KeyValValuedouble(x0p_mc_));
  bp_mc_   = keyval->doublevalue("bp_mc", KeyValValuedouble(bp_mc_));
  cp_mc_   = keyval->doublevalue("cp_mc", KeyValValuedouble(cp_mc_));
  x0f_mc_  = keyval->doublevalue("x0f_mc", KeyValValuedouble(x0f_mc_));
  bf_mc_   = keyval->doublevalue("bf_mc", KeyValValuedouble(bf_mc_));
  cf_mc_   = keyval->doublevalue("cf_mc", KeyValValuedouble(cf_mc_));
  x0_alpha_mc_ = keyval->doublevalue("x0_alpha_mc",
                                     KeyValValuedouble(x0_alpha_mc_));
  b_alpha_mc_  = keyval->doublevalue("b_alpha_mc",
                                     KeyValValuedouble(b_alpha_mc_));
  c_alpha_mc_  = keyval->doublevalue("c_alpha_mc",
                                     KeyValValuedouble(c_alpha_mc_));

  x0p_rpa_ = keyval->doublevalue("x0p_rpa", KeyValValuedouble(x0p_rpa_));
  bp_rpa_  = keyval->doublevalue("bp_rpa", KeyValValuedouble(bp_rpa_));
  cp_rpa_  = keyval->doublevalue("cp_rpa", KeyValValuedouble(cp_rpa_));
  x0f_rpa_ = keyval->doublevalue("x0f_rpa", KeyValValuedouble(x0f_rpa_));
  bf_rpa_  = keyval->doublevalue("bf_rpa", KeyValValuedouble(bf_rpa_));
  cf_rpa_  = keyval->doublevalue("cf_rpa", KeyValValuedouble(cf_rpa_));
  x0_alpha_rpa_ = keyval->doublevalue("x0_alpha_rpa",
                                      KeyValValuedouble(x0_alpha_rpa_));
  b_alpha_rpa_  = keyval->doublevalue("b_alpha_rpa",
                                      KeyValValuedouble(b_alpha_rpa_));
  c_alpha_rpa_  = keyval->doublevalue("c_alpha_rpa",
                                      KeyValValuedouble(c_alpha_rpa_));   
}

void
VWNLCFunctional::save_data_state(StateOut& s)
{
  LSDACFunctional::save_data_state(s);

  s.put(Ap_);
  s.put(Af_);
  s.put(A_alpha_);

  s.put(x0p_mc_);
  s.put(bp_mc_);
  s.put(cp_mc_);
  s.put(x0f_mc_);
  s.put(bf_mc_);
  s.put(cf_mc_);
  s.put(x0_alpha_mc_);
  s.put(b_alpha_mc_);
  s.put(c_alpha_mc_);
  
  s.put(x0p_rpa_);
  s.put(bp_rpa_);
  s.put(cp_rpa_);
  s.put(x0f_rpa_);
  s.put(bf_rpa_);
  s.put(cf_rpa_);
  s.put(x0_alpha_rpa_);
  s.put(b_alpha_rpa_);
  s.put(c_alpha_rpa_);
}

void
VWNLCFunctional::init_constants()
{
  Ap_ = 0.0310907;
  Af_ = 0.01554535;
  A_alpha_ = -1./(6.*M_PI*M_PI);

  x0p_mc_ = -0.10498;
  bp_mc_  = 3.72744;
  cp_mc_  = 12.9352;
  x0f_mc_ = -0.32500;
  bf_mc_  = 7.06042;
  cf_mc_  = 18.0578;
  x0_alpha_mc_ = -0.00475840;
  b_alpha_mc_  = 1.13107;
  c_alpha_mc_  = 13.0045;
  
  x0p_rpa_ = -0.409286;
  bp_rpa_  = 13.0720;
  cp_rpa_  = 42.7198;
  x0f_rpa_ = -0.743294;
  bf_rpa_  = 20.1231;
  cf_rpa_  = 101.578;
  x0_alpha_rpa_ = -0.228344;
  b_alpha_rpa_  = 1.06835;
  c_alpha_rpa_  = 11.4813;
}

double
VWNLCFunctional::F(double x, double A, double x0, double b, double c)
{
  double x2 = x*x;
  double x02 = x0*x0;
  double Xx = x2 + b*x + c;
  double Xx0 = x02 + b*x0 + c;
  double Q = sqrt(4.*c-b*b);
  double res
      = A * ( log(x2/Xx)
              + 2.*b/Q * atan(Q/(2.*x+b))
              - b*x0/Xx0 * ( log((x-x0)*(x-x0)/Xx)
                           + 2.*(b+2.*x0)/Q * atan(Q/(2.*x+b)) ) );
  return res;
}

double
VWNLCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
{
  double x2 = x*x;
  double x02 = x0*x0;
  double Xx = x2 + b*x +c;
  double Xx0 = x02 + b*x0 + c;
  double Q = sqrt(4.*c-b*b);
  double res
      = A * ( 1./x2 - 1./Xx - b/(2.*Xx*x) 
          + ((x0*(2.*x0+b))/Xx0 - 1) * (2.*b)/(x*(Q*Q+(2.*x+b)*(2.*x+b))) 
          - (b*x0)/(x*(x-x0)*Xx0) + (b*x0*(1+(b/(2.*x))))/(Xx0*Xx) );
  return res;
}

void
VWNLCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
}

/////////////////////////////////////////////////////////////////////////////
// VWN1LCFunctional
// Coded by Matt Leininger

static ClassDesc VWN1LCFunctional_cd(
  typeid(VWN1LCFunctional),"VWN1LCFunctional",1,"public LSDACFunctional",
  0, create<VWN1LCFunctional>, create<VWN1LCFunctional>);

VWN1LCFunctional::VWN1LCFunctional(StateIn& s):
  SavableState(s),
  VWNLCFunctional(s)
{
  s.get(x0p_);
  s.get(bp_);
  s.get(cp_);
  s.get(x0f_);
  s.get(bf_);
  s.get(cf_);
}

VWN1LCFunctional::VWN1LCFunctional()
{
  x0p_ = x0p_mc_;
  bp_  = bp_mc_;
  cp_  = cp_mc_;
  x0f_ = x0f_mc_;
  bf_  = bf_mc_;
  cf_  = cf_mc_;
}

VWN1LCFunctional::VWN1LCFunctional(int use_rpa)
{
  if (use_rpa) {
      x0p_ = x0p_rpa_;
      bp_  = bp_rpa_;
      cp_  = cp_rpa_;
      x0f_ = x0f_rpa_;
      bf_  = bf_rpa_;
      cf_  = cf_rpa_;
    }
  else {
      x0p_ = x0p_mc_;
      bp_  = bp_mc_;
      cp_  = cp_mc_;
      x0f_ = x0f_mc_;
      bf_  = bf_mc_;
      cf_  = cf_mc_;
    }
  
}

VWN1LCFunctional::VWN1LCFunctional(const Ref<KeyVal>& keyval):
  VWNLCFunctional(keyval)
{
  int vwn1rpa = keyval->booleanvalue("rpa", KeyValValueboolean(0));
  if (vwn1rpa) {
      x0p_ = x0p_rpa_;
      bp_  = bp_rpa_;
      cp_  = cp_rpa_;
      x0f_ = x0f_rpa_;
      bf_  = bf_rpa_;
      cf_  = cf_rpa_;
    }
  else {
      x0p_ = x0p_mc_;
      bp_  = bp_mc_;
      cp_  = cp_mc_;
      x0f_ = x0f_mc_;
      bf_  = bf_mc_;
      cf_  = cf_mc_;
    }

  x0p_ = keyval->doublevalue("x0p", KeyValValuedouble(x0p_));
  bp_ = keyval->doublevalue("bp", KeyValValuedouble(bp_));
  cp_ = keyval->doublevalue("cp", KeyValValuedouble(cp_));
  x0f_ = keyval->doublevalue("x0f", KeyValValuedouble(x0f_));
  bf_ = keyval->doublevalue("bf", KeyValValuedouble(bf_));
  cf_ = keyval->doublevalue("cf", KeyValValuedouble(cf_));
}

VWN1LCFunctional::~VWN1LCFunctional()
{
}

void
VWN1LCFunctional::save_data_state(StateOut& s)
{
  VWNLCFunctional::save_data_state(s);
  s.put(x0p_);
  s.put(bp_);
  s.put(cp_);
  s.put(x0f_);
  s.put(bf_);
  s.put(cf_);
}

// Based on the VWN1 functional in Vosko, Wilk, and Nusair, Can. J. Phys.
// 58, 1200, (1980).
void
VWN1LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();

  double four_thirds = 4./3.;
  double one_third = 1./3.;
  double two_thirds = 2./3.;
  double sixth = 1./6.;
  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;
  
  double epc    = F(x, Ap_, x0p_, bp_, cp_);
  double efc    = F(x, Af_, x0f_, bf_, cf_);

  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double delta_ec = f * (efc - epc);
  double ec = epc + delta_ec;

  ec_local = ec;
  od.energy = ec * rho;

  if (compute_potential_) {
    if (!spin_polarized_) {
      double depc_dr_s0 = dFdr_s(x, Ap_, x0p_, bp_, cp_);
      double dec_dr_s = depc_dr_s0;
      od.df_drho_a = od.df_drho_b = ec - (rs/3.)*dec_dr_s;
      decrs = dec_dr_s;
      deczeta = 0.;
    }
    else {
      double depc_dr_s0 = dFdr_s(x, Ap_, x0p_, bp_, cp_);
      double defc_dr_s1 = dFdr_s(x, Af_, x0f_, bf_, cf_);
      double dec_dr_s = depc_dr_s0 + f * (defc_dr_s1 - depc_dr_s0);
      double fp = two_thirds * (pow((1+zeta),one_third) 
            - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      // double dec_dzeta = depc_dr_s0 + fp * (efc - epc) + f * (defc_dr_s1 - depc_dr_s0);
      double dec_dzeta = fp * (efc - epc);       
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
      } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// VWN2LCFunctional
// Coded by Matt Leininger
static ClassDesc VWN2LCFunctional_cd(
  typeid(VWN2LCFunctional),"VWN2LCFunctional",1,"public VWNLCFunctional",
  0, create<VWN2LCFunctional>, create<VWN2LCFunctional>);

VWN2LCFunctional::VWN2LCFunctional(StateIn& s):
  SavableState(s),
  VWNLCFunctional(s)
{
}

VWN2LCFunctional::VWN2LCFunctional()
{
}

VWN2LCFunctional::VWN2LCFunctional(const Ref<KeyVal>& keyval):
  VWNLCFunctional(keyval)
{
}

VWN2LCFunctional::~VWN2LCFunctional()
{
}

void
VWN2LCFunctional::save_data_state(StateOut& s)
{
  VWNLCFunctional::save_data_state(s);
}

void
VWN2LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();
  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  const double sixth = 1./6.;
  const double four_thirds = 4./3.;
  const double one_third = 1./3.;
  const double two_thirds = 2./3.;

  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;
  
  // Monte Carlo fitting parameters 
  double epc_mc    = F(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
  double efc_mc    = F(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
  //double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
  // RPA fitting parameters
  double epc_rpa    = F(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
  double efc_rpa    = F(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
  double alphac_rpa = F(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);

  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc_rpa - epc_rpa;
  double delta_e_mc  = efc_mc - epc_mc;
  double beta = (fpp0*delta_e_rpa / alphac_rpa) - 1.;
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. + beta * zeta4);
  double delta_ec =  delta_erpa_rszeta + f*(delta_e_mc - delta_e_rpa);
  double ec = epc_mc + delta_ec;

  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      //double zeta3 = zeta2*zeta;
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, Ap_,x0p_mc_, bp_mc_, cp_mc_); 
      double defc_dr_s1_mc = dFdr_s(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
      double dalphac_dr_s_mc = dFdr_s(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
      // RPA fitting parameters
      double depc_dr_s0_rpa = dFdr_s(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
      double defc_dr_s1_rpa = dFdr_s(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
      double dalphac_dr_s_rpa = dFdr_s(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);
    
      double fp = two_thirds * (pow((1+zeta),one_third)
             - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      
      // RPA fitting parameters
      double ddelta_e_rpa = defc_dr_s1_rpa - depc_dr_s0_rpa;
      double ddeltae_rpa_dr_s = f / fpp0 * (1 - zeta4)* dalphac_dr_s_rpa 
                    + f * zeta4 * ddelta_e_rpa;
      double ddeltae_rpa_dzeta = alphac_rpa / fpp0 * 
        ( fp * (1.-zeta4) - 4.* f * zeta*zeta2) 
        + delta_e_rpa * ( fp*zeta4 + 4.*f*zeta*zeta2);
              
      // Monte Carlo fitting parameters
      double ddelta_e_mc = defc_dr_s1_rpa - depc_dr_s0_mc;
      // double dec_dr_s_mc = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
      //                   + dalphac_dr_s * f / fpp0 * (1 - zeta4);
      // double dec_dzeta_mc = 4.* zeta3 * f * (efc_mc - epc_mc - (alphac_mc/fpp0))
      //        + fp * (zeta4 * (efc_mc - epc_mc) + (1-zeta4)*(alphac_mc/fpp0));
      
      double dec_dr_s = depc_dr_s0_mc + ddeltae_rpa_dr_s + f * (ddelta_e_mc - ddelta_e_rpa);
      
      double dec_dzeta = ddeltae_rpa_dzeta + fp * (delta_e_mc - delta_e_rpa);
                       
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN3LCFunctional
// Coded by Matt Leininger
static ClassDesc VWN3LCFunctional_cd(
  typeid(VWN3LCFunctional),"VWN3LCFunctional",1,"public VWNLCFunctional",
  0, create<VWN3LCFunctional>, create<VWN3LCFunctional>);

VWN3LCFunctional::VWN3LCFunctional(StateIn& s):
  SavableState(s),
  VWNLCFunctional(s)
{
  s.get(monte_carlo_prefactor_);
  s.get(monte_carlo_e0_);
}

VWN3LCFunctional::VWN3LCFunctional(int mcp, int mce0)
{
  monte_carlo_prefactor_ = mcp;
  monte_carlo_e0_ = mce0;
}

VWN3LCFunctional::VWN3LCFunctional(const Ref<KeyVal>& keyval):
  VWNLCFunctional(keyval)
{
    monte_carlo_prefactor_ = keyval->booleanvalue("monte_carlo_prefactor",
                                                  KeyValValueboolean(1));
    monte_carlo_e0_ = keyval->booleanvalue("monte_carlo_e0",
                                           KeyValValueboolean(1));
}

VWN3LCFunctional::~VWN3LCFunctional()
{
}

void
VWN3LCFunctional::save_data_state(StateOut& s)
{
  VWNLCFunctional::save_data_state(s);
  s.put(monte_carlo_prefactor_);
  s.put(monte_carlo_e0_);
}

// based on the equations given on a NIST WWW site
void
VWN3LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();
  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  const double sixth = 1./6.;
  const double four_thirds = 4./3.;
  const double one_third = 1./3.;
  const double two_thirds = 2./3.;

  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;
  
  // Monte Carlo fitting parameters 
  double epc0_mc    = F(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
  double efc1_mc    = F(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
  // double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
  // RPA fitting parameters
  double epc0_rpa    = F(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
  double efc1_rpa    = F(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
  double alphac_rpa = F(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);
    
  double f = 9./8.*fpp0*( pow(1.+zeta, four_thirds) + pow(1.-zeta, four_thirds) - 2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc1_rpa - epc0_rpa;
  double delta_e_mc  = efc1_mc  - epc0_mc;
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1.-zeta4) + f * zeta4 * delta_e_rpa;
  double delta_ec;
  if (!monte_carlo_prefactor_) delta_ec = delta_erpa_rszeta;
  else delta_ec = delta_e_mc/delta_e_rpa * delta_erpa_rszeta;

  double ec;
  if (monte_carlo_e0_) ec = epc0_mc;
  else ec = epc0_rpa;
  ec += delta_ec;
  
  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
      double defc_dr_s1_mc = dFdr_s(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
      // double dalphac_dr_s_mc = dFdr_s(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
      // RPA fitting parameters
      double depc_dr_s0_rpa = dFdr_s(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
      double defc_dr_s1_rpa = dFdr_s(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
      double dalphac_dr_s_rpa = dFdr_s(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);

      double fp = two_thirds * (pow((1+zeta),one_third)
                  - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
           
      // RPA fitting parameters
      double ddelta_e_rpa = defc_dr_s1_rpa - depc_dr_s0_rpa;
      double ddeltae_rpa_dr_s = f / fpp0 * (1 - zeta4)* dalphac_dr_s_rpa 
                    + f * zeta4 * ddelta_e_rpa;
      double ddeltae_rpa_dzeta = alphac_rpa / fpp0 * 
         ( fp * (1.-zeta4) - 4.* f * zeta*zeta2) 
      + delta_e_rpa * ( fp*zeta4 + 4.*f*zeta*zeta2 );

      // Monte Carlo fitting parameters
      double ddelta_e_mc = defc_dr_s1_mc - depc_dr_s0_mc;
      
      double dec_dzeta, dec_dr_s;
      if (!monte_carlo_prefactor_) {
        dec_dzeta = ddeltae_rpa_dzeta;
        if (monte_carlo_e0_) dec_dr_s = depc_dr_s0_mc;
        else dec_dr_s = depc_dr_s0_rpa;
        dec_dr_s += ddeltae_rpa_dr_s;
        }
      else {
        dec_dzeta = delta_e_mc / delta_e_rpa * ddeltae_rpa_dzeta;
        if (monte_carlo_e0_) dec_dr_s = depc_dr_s0_mc;
        else dec_dr_s = depc_dr_s0_rpa;
        dec_dr_s += delta_erpa_rszeta/delta_e_rpa * ddelta_e_mc
                   + delta_e_mc/delta_e_rpa * ddeltae_rpa_dr_s 
                   - delta_erpa_rszeta*delta_e_mc/(delta_e_rpa*delta_e_rpa)
                   * ddelta_e_rpa;
        }        
        
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN4LCFunctional
// Coded by Matt Leininger
static ClassDesc VWN4LCFunctional_cd(
  typeid(VWN4LCFunctional),"VWN4LCFunctional",1,"public VWNLCFunctional",
  0, create<VWN4LCFunctional>, create<VWN4LCFunctional>);

VWN4LCFunctional::VWN4LCFunctional(StateIn& s):
  SavableState(s),
  VWNLCFunctional(s)
{
  s.get(monte_carlo_prefactor_);
}

VWN4LCFunctional::VWN4LCFunctional()
{
  monte_carlo_prefactor_ = 0;
}

VWN4LCFunctional::VWN4LCFunctional(const Ref<KeyVal>& keyval):
  VWNLCFunctional(keyval)
{
  monte_carlo_prefactor_ = keyval->booleanvalue("monte_carlo_prefactor",
                                                KeyValValueboolean(0));
}

VWN4LCFunctional::~VWN4LCFunctional()
{
}

void
VWN4LCFunctional::save_data_state(StateOut& s)
{
  VWNLCFunctional::save_data_state(s);
  s.put(monte_carlo_prefactor_);
}

// based on the equations given on a NIST WWW site
void
VWN4LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();
  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  const double sixth = 1./6.;
  const double four_thirds = 4./3.;
  const double one_third = 1./3.;
  const double two_thirds = 2./3.;

  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;
  
  // Monte Carlo fitting parameters 
  double epc_mc    = F(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
  double efc_mc    = F(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
  // double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
  // RPA fitting parameters
  double epc_rpa    = F(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
  double efc_rpa    = F(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
  double alphac_rpa = F(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);

  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc_rpa - epc_rpa;
  double delta_e_mc  = efc_mc - epc_mc;
  // double beta = fpp0 * delta_e_mc / alphac_rpa - 1.;
  // double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. + beta * zeta4);
  // use delta_e_rpa here instead of delta_e_mc to get VWN3
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. - zeta4)
                           + f * delta_e_mc * zeta4;
  double delta_ec;
  if (!monte_carlo_prefactor_) delta_ec = delta_erpa_rszeta;
  else delta_ec = delta_e_mc/delta_e_rpa * delta_erpa_rszeta;
  // double ec = epc_rpa + delta_ec;
  double ec = epc_mc + delta_ec;
  
  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      double zeta3 = zeta2*zeta;
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
      double defc_dr_s1_mc = dFdr_s(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
      // double dalphac_dr_s_mc = dFdr_s(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
      // RPA fitting parameters
      double depc_dr_s0_rpa = dFdr_s(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
      double defc_dr_s1_rpa = dFdr_s(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
      double dalphac_dr_s_rpa = dFdr_s(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_,
                                       c_alpha_rpa_);
    
      double fp = two_thirds * (pow((1+zeta),one_third)
             - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      
      double ddelta_e_rpa = defc_dr_s1_rpa - depc_dr_s0_rpa;
      double ddelta_e_mc = defc_dr_s1_mc - depc_dr_s0_mc;
      // use ddelta_e_rpa here instead of ddelta_e_mc to get VWN3
      double ddeltae_rpa_dr_s = f / fpp0 * (1 - zeta4)* dalphac_dr_s_rpa 
                    + f * zeta4 * ddelta_e_mc;
      // use ddelta_e_rpa here instead of ddelta_e_mc to get VWN3
      double ddeltae_rpa_dzeta = alphac_rpa / fpp0 * 
        ( fp * (1.-zeta4) - 4.* f * zeta3) 
        + delta_e_mc * ( fp*zeta4 + 4.*f*zeta3);

      double dec_dzeta, dec_dr_s;
      if (!monte_carlo_prefactor_) {
          // dec_dr_s = depc_dr_s0_rpa + ddeltae_rpa_dr_s;
          dec_dr_s = depc_dr_s0_mc + ddeltae_rpa_dr_s;
          dec_dzeta = ddeltae_rpa_dzeta;
        }
      else {
          dec_dr_s = depc_dr_s0_rpa + delta_erpa_rszeta/delta_e_rpa * ddelta_e_mc
                   + delta_e_mc/delta_e_rpa * ddeltae_rpa_dr_s 
                   - delta_erpa_rszeta*delta_e_mc/(delta_e_rpa*delta_e_rpa)* ddelta_e_rpa;
          dec_dzeta = delta_e_mc / delta_e_rpa * ddeltae_rpa_dzeta;      
        }
      
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1.)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1.)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN5LCFunctional
// Coded by Matt Leininger

static ClassDesc VWN5LCFunctional_cd(
  typeid(VWN5LCFunctional),"VWN5LCFunctional",1,"public VWNLCFunctional",
  0, create<VWN5LCFunctional>, create<VWN5LCFunctional>);

VWN5LCFunctional::VWN5LCFunctional(StateIn& s):
  SavableState(s),
  VWNLCFunctional(s)
{
}

VWN5LCFunctional::VWN5LCFunctional()
{
}

VWN5LCFunctional::VWN5LCFunctional(const Ref<KeyVal>& keyval):
  VWNLCFunctional(keyval)
{
}

VWN5LCFunctional::~VWN5LCFunctional()
{
}

void
VWN5LCFunctional::save_data_state(StateOut& s)
{
  VWNLCFunctional::save_data_state(s);
}

// based on the equations given on a NIST WWW site
// based on the VWN5 functional in Vosko, Wilk, and Nusair, Can. J. Phys.
// 58, 1200, (1980). and Perdew and Wang, Phys. Rev. B, 45, 13244, (1992).
void
VWN5LCFunctional::point_lc(const PointInputData &id, PointOutputData &od, 
                          double &ec_local, double &decrs, double &deczeta)
{
  od.zero();
  const double fpp0 = 4./9. * 1./(pow(2., (1./3.)) - 1.);
  const double sixth = 1./6.;
  const double four_thirds = 4./3.;
  const double one_third = 1./3.;
  const double two_thirds = 2./3.;

  double rho = id.a.rho + id.b.rho;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double x = pow(3./(4.*M_PI*rho), sixth);
  double rs = x*x;

  double epc    = F(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
  double efc    = F(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
  double alphac = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
     
  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double beta = fpp0 * (efc - epc) / alphac - 1.;
  double delta_ec = alphac * f / fpp0 * (1. + beta * zeta4);
  double ec = epc + delta_ec;

  ec_local = ec;
  od.energy = ec * rho;

  if (compute_potential_) {
    if (!spin_polarized_) {
      double depc_dr_s0 = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_); 
      double dec_dr_s = depc_dr_s0;
      od.df_drho_a = od.df_drho_b = ec - (rs/3.)*dec_dr_s;
      decrs = dec_dr_s;
      deczeta = 0.;
    }
    else {
      double zeta3 = zeta2*zeta;
      double depc_dr_s0 = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
      double defc_dr_s1 = dFdr_s(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
      double dalphac_dr_s = dFdr_s(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
      double dec_dr_s = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
                        + dalphac_dr_s * f / fpp0 * (1 - zeta4);
      double fp = two_thirds * (pow((1+zeta),one_third) 
            - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      double dec_dzeta = 4.* zeta3 * f * (efc - epc - (alphac/fpp0))
              + fp * (zeta4 * (efc - epc) + (1-zeta4)*(alphac/fpp0));
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
      } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// XalphaFunctional

static ClassDesc XalphaFunctional_cd(
  typeid(XalphaFunctional),"XalphaFunctional",1,"public DenFunctional",
  0, create<XalphaFunctional>, create<XalphaFunctional>);

XalphaFunctional::XalphaFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(alpha_);
  factor_ = alpha_ * 2.25 * pow(3.0/(4.*M_PI), 1.0/3.0);
}

XalphaFunctional::XalphaFunctional()
{
  alpha_ = 0.70;
  factor_ = alpha_ * 2.25 * pow(3.0/(4.*M_PI), 1.0/3.0);
}

XalphaFunctional::XalphaFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  alpha_ = keyval->doublevalue("alpha", KeyValValuedouble(0.70));
  factor_ = alpha_ * 2.25 * pow(3.0/(4.*M_PI), 1.0/3.0);
}

XalphaFunctional::~XalphaFunctional()
{
}

void
XalphaFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(alpha_);
}

void
XalphaFunctional::point(const PointInputData &id,
                        PointOutputData &od)
{
  od.zero();
  const double four_thirds = 4./3.;
  
  if (!spin_polarized_) {
      od.energy = - 2.0 * factor_ * id.a.rho * id.a.rho_13;
      if (compute_potential_) {
          od.df_drho_a = -four_thirds * factor_ * id.a.rho_13;
          od.df_drho_b = od.df_drho_a;
              }
    }
  else {
      double rhoa43 = id.a.rho * id.a.rho_13;
      double rhob43 = id.b.rho * id.b.rho_13;
      od.energy = - factor_ * (rhoa43 + rhob43);
      if (compute_potential_) {
          od.df_drho_a = -four_thirds * factor_ * id.a.rho_13;
          od.df_drho_b = -four_thirds * factor_ * id.b.rho_13;
              }
    }
}

void
XalphaFunctional::print(std::ostream& o) const
{
  o
    << indent << scprintf("XalphaFunctional: alpha = %12.8f", alpha_) << std::endl;
}

/////////////////////////////////////////////////////////////////////////////
// Becke88XFunctional

static ClassDesc Becke88XFunctional_cd(
  typeid(Becke88XFunctional),"Becke88XFunctional",1,"public DenFunctional",
  0, create<Becke88XFunctional>, create<Becke88XFunctional>);

Becke88XFunctional::Becke88XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(beta_);
  beta6_ = 6. * beta_;
}

Becke88XFunctional::Becke88XFunctional()
{
  beta_ = 0.0042;
  beta6_ = 6. * beta_;
}

Becke88XFunctional::Becke88XFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  beta_ = keyval->doublevalue("beta", KeyValValuedouble(0.0042));
  beta6_ = 6. * beta_;
}

Becke88XFunctional::~Becke88XFunctional()
{
}

void
Becke88XFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(beta_);
}

int
Becke88XFunctional::need_density_gradient()
{
  return 1;
}

// Becke's exchange
// From:  C.W. Murray et al.  Mol. Phys. Vol 78  pp 997-1014 (1993)
// originally coded by Mike Colvin
void
Becke88XFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();
  // Preset terms from Murray's paper
  const double beta = beta_;
  const double beta6 = beta6_;

  // const double beta6=0.0252;

  // Use simplified formula
  double rho_a_13 = pow(id.a.rho,(1./3.));
  double rho_a_43 = id.a.rho*rho_a_13;
  double xa, xa2, ga_denom, ga_denom2;
  if (id.a.rho > MIN_DENSITY) {
    xa = sqrt(id.a.gamma)/rho_a_43;
    xa2 = xa*xa;
    ga_denom = 1/(1.+beta6*xa*asinh(xa));
    ga_denom2 = ga_denom*ga_denom;
    }
  else {
    xa = xa2 = 0.;
    ga_denom = ga_denom2 = 1.;
    }
  double Fa = sqrt(1.+xa2);
  double Ha = 1. - 6.*beta*xa2/Fa;
  double ex;
  if (id.a.rho > MIN_DENSITY) ex = -rho_a_43*beta*xa2*ga_denom;
  else ex = 0.;
 
  if (compute_potential_) {
      if (id.a.rho > MIN_DENSITY) {
        od.df_drho_a = 4./3. * beta * rho_a_13 * xa2 * ga_denom2 * Ha;
        od.df_dgamma_aa = -beta * ga_denom / (2.*rho_a_43) * (1. + ga_denom*Ha);
        }
      else od.df_drho_a = od.df_dgamma_aa = 0.;
      od.df_drho_b=od.df_drho_a;
      od.df_dgamma_bb=od.df_dgamma_aa;
         }

  if (spin_polarized_) {
      double rho_b_13 = pow(id.b.rho,(1./3.));
      double rho_b_43 = id.b.rho*rho_b_13;
      double xb, xb2, gb_denom, gb_denom2;
      if (id.b.rho > MIN_DENSITY) {
        xb = sqrt(id.b.gamma)/rho_b_43;
        xb2 = xb*xb;
        gb_denom = 1./(1.+beta6*xb*asinh(xb));
        gb_denom2 = gb_denom*gb_denom;
        }
      else {
        xb = xb2 = 0.;
        gb_denom = gb_denom2 = 1.;
        }
      double Fb = sqrt(1.+xb2);
      double Hb = 1. - 6.*beta*xb2/Fb;   
      ex += -rho_b_43*beta*xb2*gb_denom;

      if (compute_potential_) {
        if (id.b.rho > MIN_DENSITY) {
          od.df_drho_b = 4./3. * beta * rho_b_13 * xb2 * gb_denom2 * Hb;
          od.df_dgamma_bb = -beta / (2.*rho_b_43) * (gb_denom + gb_denom2*Hb);        
          }
        else od.df_drho_b = od.df_dgamma_bb = 0.;
      }
    }
  else ex += ex;

  od.energy = ex;
}


/////////////////////////////////////////////////////////////////////////////
// LYPCFunctional
// Coded by Matt Leininger

static ClassDesc LYPCFunctional_cd(
  typeid(LYPCFunctional),"LYPCFunctional",1,"public DenFunctional",
  0, create<LYPCFunctional>, create<LYPCFunctional>);

LYPCFunctional::LYPCFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(a_);
  s.get(b_);
  s.get(c_);
  s.get(d_);
}

LYPCFunctional::LYPCFunctional()
{
  init_constants();
}

LYPCFunctional::LYPCFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  a_ = keyval->doublevalue("a", KeyValValuedouble(a_));
  b_ = keyval->doublevalue("b", KeyValValuedouble(b_));
  c_ = keyval->doublevalue("c", KeyValValuedouble(c_));
  d_ = keyval->doublevalue("d", KeyValValuedouble(d_));
}

LYPCFunctional::~LYPCFunctional()
{
}

void
LYPCFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(a_);
  s.put(b_);
  s.put(c_);
  s.put(d_);
}

void
LYPCFunctional::init_constants()
{
  a_ = 0.04918;
  b_ = 0.132;
  c_ = 0.2533;
  d_ = 0.349;
}

int
LYPCFunctional::need_density_gradient()
{
  return 1;
}

// Lee-Yang-Parr correlation
// From: Burkhard Miehlich, et al.  Chem Phys. Lett. Vol. 157 pp200-206 (1989)
// Original LYP paper Phys. Rev. B Vol. 37 pp785-789 (1988)
// originally coded by Mike Colvin
// potential terms added by Matt Leininger
void
LYPCFunctional::point(const PointInputData &id,
                     PointOutputData &od)
{
  od.zero();
  double ec;

  // Precalculate terms for efficiency
  double dens=id.a.rho+id.b.rho;
  double dens2=dens*dens;
  // double grad=id.a.gamma+id.b.gamma;
  // double grad=sqrt( pow(id.a.del_rho[0]+id.b.del_rho[0],2.) 
  //                + pow(id.a.del_rho[1]+id.b.del_rho[1],2.) 
  //                  + pow(id.a.del_rho[2]+id.b.del_rho[2],2.) );
  // double grad2=grad*grad;
  double dens1_3=pow(dens,-1./3.);

  // Precalculate terms defined in Miehlich's paper
  const double a = a_;
  const double b = b_;
  const double c = c_;
  const double d = d_;
  double omega=exp(-c*dens1_3)/(1.+d*dens1_3)*pow(dens,-11./3.);
  double delta=c*dens1_3+d*dens1_3/(1.+d*dens1_3);
  double cf=0.3*pow(3.* M_PI*M_PI,2./3.);
  double denom = 1.+d*dens1_3;
  
  double dens_a2=id.a.rho*id.a.rho;
  double dens_b2=id.b.rho*id.b.rho;
  double dens_ab=id.a.rho*id.b.rho;
  double grad_a2=id.a.gamma;
  double grad_b2=id.b.gamma;
  double grad_ab=id.gamma_ab;

  double eflyp_1 = -4.*a*dens_ab/(dens*denom);
  double intermediate_1 = pow(2.,2./3.)*144.*cf*(pow(id.a.rho,8./3.)+pow(id.b.rho,8./3.))
           + (47.-7.*delta)*(grad_a2+grad_b2+2.*grad_ab)-(45.-delta)*(grad_a2+grad_b2)
           + 2.*(11.-delta)/dens*(id.a.rho*grad_a2+id.b.rho*grad_b2); 
  double intermediate_2 = -4./3.*dens2*grad_ab - (dens_a2*grad_b2+dens_b2*grad_a2);
  double intermediate_3 = dens_ab/18.* intermediate_1 + intermediate_2;  
  double eflyp_2 = -omega*a*b*intermediate_3; 
  ec = eflyp_1 + eflyp_2;
  od.energy = ec;
 
    if (compute_potential_) {
       double dens4_3 = pow(dens,-4./3.);
       double expo = exp(-c*dens1_3);
       double ddelta_drho_a = 1./3* (d*d*dens4_3*dens1_3/(denom*denom)
                                     - delta/dens);
       double domega_drho_a = -1./3.*omega*dens4_3*(11./dens1_3 - c - d/denom);
       double df1_drho_a;
       df1_drho_a = -4.*a*id.b.rho/(dens*denom) *
                    (id.a.rho/3.*d*dens4_3/denom + 1. - id.a.rho/dens);
       // if (id.a.rho > MIN_DENSITY) 
       //  df1_drho_a = -4.*a*dens_ab/(dens*denom)
       //                  * (1./3.*d*dens4_3/denom + 1/id.a.rho - 1./dens); 
       // else df1_drho_a = 0.;
       double df2_drho_a = -domega_drho_a*a*b*intermediate_3 
                  - omega*a*b*( id.b.rho/18.* intermediate_1
                  + dens_ab/18.*(144.*pow(2.,2./3.)*cf*8./3.*pow(id.a.rho,5./3.)
                  +  2.*(11.-delta)*grad_a2/dens 
                  - 2./dens*ddelta_drho_a*(id.a.rho*grad_a2 + id.b.rho*grad_b2) 
                  - 2.*(11.-delta)/dens2*(id.a.rho*grad_a2 + id.b.rho*grad_b2) 
                  - 7.*ddelta_drho_a*(grad_a2+grad_b2+2.*grad_ab)
                  + ddelta_drho_a*(grad_a2+grad_b2) )  
                  - 8./3.*dens*grad_ab - 2.*id.a.rho*grad_b2 );
        
       od.df_drho_a = df1_drho_a + df2_drho_a;
       
       od.df_dgamma_aa = -omega*a*b  
                * (dens_ab/9.*(1.-3.*delta + id.a.rho*(11.-delta)/dens) - dens_b2);
       od.df_dgamma_ab = -omega*a*b*(dens_ab/9.*(47.-7.*delta) - 4./3.*dens2);
       od.df_drho_b = od.df_drho_a;
       od.df_dgamma_bb = od.df_dgamma_aa;
       if (spin_polarized_) {
          double ddelta_drho_b = ddelta_drho_a;
          double domega_drho_b = domega_drho_a;
          double df1_drho_b;
          df1_drho_b = -4.*a*id.a.rho/(dens*denom) *
                       (id.b.rho/3.*d*dens4_3/denom + 1. - id.b.rho/dens);
          //if (id.b.rho > MIN_DENSITY)
          //  df1_drho_b = -4.*a*dens_ab/(dens*denom)
          //                  * (1./3.*d*dens4_3/denom + 1./id.b.rho - 1./dens);
          //else df1_drho_b = 0.;
          double df2_drho_b = -domega_drho_b*a*b*intermediate_3 
                  - omega*a*b*( id.a.rho/18.* intermediate_1 
                  + dens_ab/18.*(144.*pow(2.,2./3.)*cf*8./3.*pow(id.b.rho,5./3.)
                  +  2.*(11.-delta)*grad_b2/dens
                  - 2./dens*ddelta_drho_b*(id.a.rho*grad_a2 + id.b.rho*grad_b2)
                  - 2.*(11.-delta)/dens2*(id.a.rho*grad_a2 + id.b.rho*grad_b2)
                  - 7.*ddelta_drho_b*(grad_a2+grad_b2+2.*grad_ab)
                  + ddelta_drho_b*(grad_a2+grad_b2) )
                  - 8./3.*dens*grad_ab - 2.*id.b.rho*grad_a2 );
          od.df_drho_b = df1_drho_b + df2_drho_b;
          od.df_dgamma_bb = -omega*a*b
                   * (dens_ab/9.*(1.-3.*delta + id.b.rho*(11.-delta)/dens) - dens_a2);
       }
    }

}

/////////////////////////////////////////////////////////////////////////////
// Perdew 1986 (P86) Correlation Functional
// J. P. Perdew, PRB, 33, 8822, 1986.
// C. W. Murray, N. C. Handy and G. J. Laming, Mol. Phys., 78, 997, 1993.
// 
// Coded by Matt Leininger

static ClassDesc P86CFunctional_cd(
  typeid(P86CFunctional),"P86CFunctional",1,"public DenFunctional",
  0, create<P86CFunctional>, create<P86CFunctional>);

P86CFunctional::P86CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(a_);
  s.get(C1_);
  s.get(C2_);
  s.get(C3_);
  s.get(C4_);
  s.get(C5_);
  s.get(C6_);
  s.get(C7_);
}

P86CFunctional::P86CFunctional()
{
  init_constants();
}

P86CFunctional::P86CFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  a_ = keyval->doublevalue("a", KeyValValuedouble(a_));
  C1_ = keyval->doublevalue("C1", KeyValValuedouble(C1_));
  C2_ = keyval->doublevalue("C2", KeyValValuedouble(C2_));
  C3_ = keyval->doublevalue("C3", KeyValValuedouble(C3_));
  C4_ = keyval->doublevalue("C4", KeyValValuedouble(C4_));
  C5_ = keyval->doublevalue("C5", KeyValValuedouble(C5_));
  C6_ = keyval->doublevalue("C6", KeyValValuedouble(C6_));
  C7_ = keyval->doublevalue("C7", KeyValValuedouble(C7_));
}

P86CFunctional::~P86CFunctional()
{
}

void
P86CFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(a_);
  s.put(C1_);
  s.put(C2_);
  s.put(C3_);
  s.put(C4_);
  s.put(C5_);
  s.put(C6_);
  s.put(C7_);
}

void
P86CFunctional::init_constants()
{
  a_ = 1.745*0.11;
  C1_ = 0.001667;
  C2_ = 0.002568;
  C3_ = 0.023266;
  C4_ = 7.389e-6;
  C5_ = 8.723;
  C6_ = 0.472;
  C7_ = 1e4*C4_;
}

int
P86CFunctional::need_density_gradient()
{
  return 1;
}

void
P86CFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();
  // Precalculate terms for efficiency
  double rho = id.a.rho + id.b.rho;
  double rho76 = pow(rho, (7./6.));
  double rho43 = pow(rho, (4./3.));
  double rho13 = pow(rho, (1./3.));
  double rho16 = pow(rho, (1./6.));
  double rs = pow( (3./(4.*M_PI*rho)), (1./3.));
  double rs2 = rs*rs;
  double rs3 = rs2*rs;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double fzeta = pow(2.,(1./3.)) *
                sqrt( pow(((1.+zeta)/2.),(5./3.)) + pow(((1.-zeta)/2.),(5./3.)) );
  double C_infin = C1_ + C2_;
  double numer = C2_ + C3_*rs + C4_*rs2;
  double denom = 1.+C5_*rs+C6_*rs2+C7_*rs3;
  double C_rho = C1_ + numer/denom;
  double gamma_aa = id.a.gamma;
  double gamma_bb = id.b.gamma;
  double gamma_ab = id.gamma_ab;
  double gamma_total = sqrt(gamma_aa + gamma_bb + 2.*gamma_ab);
  double gamma_total2 = gamma_total*gamma_total;
  double Phi = a_ * C_infin * gamma_total / (C_rho * rho76);
  double fp86 = exp(-Phi) * C_rho * gamma_total2 / (fzeta*rho43);
  if (rho < MIN_DENSITY) fp86 = 0.;    
  od.energy = fp86;

  if (compute_potential_) {
      double drs_drhoa = -rs/(3.*rho);
      double dCrho_drhoa = drs_drhoa/denom *
                           (C3_+2.*C4_*rs - numer/denom * (C5_+2.*C6_*rs+3.*C7_*rs2));
      double dCrho_drhob = dCrho_drhoa;
      double dzeta_drhoa = 1./rho * (1.-zeta);
      double dzeta_drhob = 1./rho * (-1.-zeta);
      double dPhi_drhoa = -Phi/(C_rho*rho76)*(dCrho_drhoa*rho76 + C_rho*(7./6.)*rho16);
      double dPhi_drhob = dPhi_drhoa;
      double dfzeta_drhoa = pow(2., (-1./3.))*1./fzeta *
      (5./3. * pow(((1.+zeta)/2.), (2./3.))*0.5*dzeta_drhoa
       + 5./3.*pow(((1.-zeta)/2.), (2./3.))*-0.5*dzeta_drhoa);
      double dfzeta_drhob = pow(2., (-1./3.))*1./fzeta *
      (5./3. * pow(((1.+zeta)/2.), (2./3.))*0.5*dzeta_drhob
       + 5./3.*pow(((1.-zeta)/2.), (2./3.))*-0.5*dzeta_drhob);
      double dfp86_drhoa = fp86/C_rho*(-dPhi_drhoa*C_rho + dCrho_drhoa)
      - fp86/(fzeta*rho43)*(dfzeta_drhoa*rho43 + fzeta*(4./3.)*rho13);
      double dfp86_drhob = fp86/C_rho*(-dPhi_drhob*C_rho + dCrho_drhob)
      - fp86/(fzeta*rho43)*(dfzeta_drhob*rho43 + fzeta*(4./3.)*rho13);
      if (rho < MIN_DENSITY) dfp86_drhoa = dfp86_drhob = 0.;
      od.df_drho_a = dfp86_drhoa;
      od.df_drho_b = dfp86_drhob;

      // gamma part of potential
      // double dPhi_dgamma_aa = Phi/(2.*gamma_total2);
      // double dPhi_dgamma_bb = dPhi_dgamma_aa;
      // double dfp86_dgamma_aa = fp86*(1./gamma_total2 - dPhi_dgamma_aa);
      double prefactor = exp(-Phi)*C_rho/( fzeta*rho43 );
      double dfp86_dgamma_aa = prefactor * (1.-Phi/2.);
      double dfp86_dgamma_bb = dfp86_dgamma_aa;
      if (rho < MIN_DENSITY) dfp86_dgamma_aa = dfp86_dgamma_bb = 0.;
      od.df_dgamma_aa = dfp86_dgamma_aa;
      od.df_dgamma_bb = dfp86_dgamma_bb;

      // double dPhi_dgamma_ab = 2.*dPhi_dgamma_aa;
      // double dfp86_dgamma_ab = fp86*(2./gamma_total2 - dPhi_dgamma_ab);
      double dfp86_dgamma_ab = prefactor * (2.-Phi);
      if (rho < MIN_DENSITY) dfp86_dgamma_ab= 0.;
      od.df_dgamma_ab = dfp86_dgamma_ab;

   }   
}


/////////////////////////////////////////////////////////////////////////////
// Perdew 1986 (P86) Correlation Functional
// J. P. Perdew, PRB, 33, 8822, 1986.
// C. W. Murray, N. C. Handy and G. J. Laming, Mol. Phys., 78, 997, 1993.
// 
// Coded by Matt Leininger

static ClassDesc NewP86CFunctional_cd(
  typeid(NewP86CFunctional),"NewP86CFunctional",1,"public DenFunctional",
  0, create<NewP86CFunctional>, create<NewP86CFunctional>);

NewP86CFunctional::NewP86CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(a_);
  s.get(C1_);
  s.get(C2_);
  s.get(C3_);
  s.get(C4_);
  s.get(C5_);
  s.get(C6_);
  s.get(C7_);
}

NewP86CFunctional::NewP86CFunctional()
{
  init_constants();
}

NewP86CFunctional::NewP86CFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  a_ = keyval->doublevalue("a", KeyValValuedouble(a_));
  C1_ = keyval->doublevalue("C1", KeyValValuedouble(C1_));
  C2_ = keyval->doublevalue("C2", KeyValValuedouble(C2_));
  C3_ = keyval->doublevalue("C3", KeyValValuedouble(C3_));
  C4_ = keyval->doublevalue("C4", KeyValValuedouble(C4_));
  C5_ = keyval->doublevalue("C5", KeyValValuedouble(C5_));
  C6_ = keyval->doublevalue("C6", KeyValValuedouble(C6_));
  C7_ = keyval->doublevalue("C7", KeyValValuedouble(C7_));
}

NewP86CFunctional::~NewP86CFunctional()
{
}

void
NewP86CFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(a_);
  s.put(C1_);
  s.put(C2_);
  s.put(C3_);
  s.put(C4_);
  s.put(C5_);
  s.put(C6_);
  s.put(C7_);
}

void
NewP86CFunctional::init_constants()
{
  a_ = 1.745*0.11;
  C1_ = 0.001667;
  C2_ = 0.002568;
  C3_ = 0.023266;
  C4_ = 7.389e-6;
  C5_ = 8.723;
  C6_ = 0.472;
  C7_ = 1e4*C4_;
}

int
NewP86CFunctional::need_density_gradient()
{
  return 1;
}

double
NewP86CFunctional::rho_deriv(double rho_a, double rho_b, double mdr)
{
  double ra0,ra1,ra2,ra3,ra4,ra5,ra6,ra7,ra8,ra9,ra10,ra11,ra12,ra13,ra14,
      ra15,ra16,ra17,ra18,ra19,ra20;
  double c_infin, c_1, c_2, c_3, c_4, c_5, c_6, c_7;
  
  c_infin = C1_ + C2_;
  c_1 = C1_; c_2 = C2_; c_3 = C3_; c_4 = C4_; c_5 = C5_; c_6 = C6_; c_7 = C7_;
 
  ra0 = pow(mdr,2.0);
  ra1 = rho_b+rho_a;
  ra2 = 1/pow(ra1,1.3333333333333333);
  ra3 = -(1.0*rho_b)+rho_a;
  ra4 = 1/ra1;
//  ra4 = 1/pow(ra1,1.0);
  ra5 = -(1.0*ra3*ra4)+1.0;
  ra6 = ra3*ra4+1.0;
  ra7 = 0.3149802624737183*pow(ra6,1.6666666666666667)+
        0.3149802624737183*pow(ra5,1.6666666666666667);
  ra8 = 1/pow(ra7,0.5);
  ra9 = 1/(ra1*ra1);
  ra10 = 1/pow(ra1,1.6666666666666667);
  ra11 = 1/pow(ra1,0.66666666666666663);
  ra12 = 1/pow(ra1,0.33333333333333331);
  ra13 = 0.62035049089940009*c_3*ra12+0.38483473155912662*c_4*ra11+c_2;
  ra14 = 0.62035049089940009*c_5*ra12+0.38483473155912662*c_6*ra11
        +0.238732414637843*c_7*ra4+1.0;
  ra15 = 1/ra14;
  ra16 = (-(0.20678349696646667*c_3*ra2)-(0.25655648770608441*c_4*ra10))*ra15-
         ((1.0*(-(0.20678349696646667*c_5*ra2)-(0.25655648770608441*c_6*ra10)-
                (0.238732414637843*c_7*ra9))*ra13)/pow(ra14,2.0));
  ra17 = 1/pow(ra1,1.1666666666666667);
  ra18 = ra13*ra15+c_1;
  ra19 = 1/ra18;
  ra20 = exp(-1.0*a_*c_infin*mdr*ra17*ra19);
  double dp86c_drho_a = 0.79370052598409979*ra0*ra2*ra8*ra18*
                        ((1.1666666666666667*a_*c_infin*mdr*ra19)/pow(ra1,2.1666666666666665)
                         +(a_*c_infin*mdr*ra17*ra16)/ra18*ra18)*
                        ra20-((1.0582673679787997*ra0*ra8*ra18*ra20)/
                              pow(ra1,2.3333333333333335))-
        ((0.3968502629920499*ra0*ra2*(0.52496710412286385*(ra4-(1.0*ra3*ra9))*
        pow(ra6,0.66666666666666663)+0.52496710412286385*(-(1.0*ra4)+ra3*ra9)
       *pow(ra5,0.66666666666666663))*ra18*ra20)/pow(ra7,1.5))+
                        0.79370052598409979*ra0*ra2*ra8*ra16*ra20;

  return dp86c_drho_a;
}

double
NewP86CFunctional::gab_deriv(double rho_a, double rho_b, double mdr)
{
  double c_infin, c_1, c_2, c_3, c_4, c_5, c_6, c_7;
  
  c_infin = C1_ + C2_;
  c_1 = C1_; c_2 = C2_; c_3 = C3_; c_4 = C4_; c_5 = C5_; c_6 = C6_; c_7 = C7_;

  double g0,g1,g2,g3,g4,g5,g6,g7;
  g0 = rho_b+rho_a;
  g1 = -(1.0*rho_b)+rho_a;
  g2 = 1/g0;
  g3 = 1/pow(0.3149802624737183*pow(g1*g2+1.0,1.6666666666666667
      )+0.3149802624737183*pow(-(1.0*g1*g2)+1.0,
                               1.6666666666666667),0.5);
  g4 = 1/pow(g0,0.66666666666666663);
  g5 = 1/pow(g0,0.33333333333333331);
  g6 = (0.62035049089940009*c_3*g5+0.38483473155912662*c_4*g4+c_2)
      /pow(0.62035049089940009*c_5*g5+0.38483473155912662*c_6*g4+
           0.238732414637843*c_7*g2+1.0,1.0)+c_1;
  g7 = exp(-1.0*a_*c_infin*mdr*1.0/pow(g0,1.1666666666666667)*1.0/g6);
  double dp86c_dmdr = (1.5874010519681996*mdr*g3*g6*g7)/pow(g0,1.3333333333333333)
                     -((0.79370052598409979*a_*c_infin*pow(mdr,2.0)*g3*g7)/pow(g0,2.5));

  return dp86c_dmdr/mdr;
}

void
NewP86CFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  // Precalculate terms for efficiency
  double rho = id.a.rho + id.b.rho;
  double rs = pow( (3./(4.*M_PI*rho)), (1./3.));
  double zet = (id.a.rho - id.b.rho)/rho;
  double c_infin = C1_ + C2_;
  double gamma_aa = id.a.gamma;
  double gamma_bb = id.b.gamma;
  double gamma_ab = id.gamma_ab;
  double mdr = sqrt(gamma_aa + gamma_bb + 2.*gamma_ab);

  double e0 = 1/pow(rho,0.66666666666666663);
  double e1 = 1/pow(rho,0.33333333333333331);
  double e2 = (0.62035049089940009*C3_*e1+0.38483473155912662*C4_*e0+C2_)
     /(0.62035049089940009*C5_*e1+0.38483473155912662*C6_*e0+(
         0.238732414637843*C7_)/rho+1.0)+C1_;
  double fp86 = (0.79370052598409979*pow(mdr,2.0)*e2)*exp(-1.0*a_*
     c_infin*mdr*1.0/e2*1.0/pow(rho,1.1666666666666667)
     )/pow(rho,1.3333333333333333)/pow(0.3149802624737183*pow(zet+
     1.0,1.6666666666666667)+0.3149802624737183*pow(-(1.0*zet)+
     1.0,1.6666666666666667),0.5);
  od.energy += fp86;

  if (compute_potential_) {
      double dfp86_drhoa = rho_deriv(id.a.rho, id.b.rho, mdr);
      double dfp86_drhob = rho_deriv(id.b.rho, id.a.rho, mdr);
      double dfp86_dgab = gab_deriv(id.a.rho, id.b.rho, mdr);
      
      od.df_drho_a += dfp86_drhoa;
      od.df_drho_b += dfp86_drhob;
      od.df_dgamma_ab += dfp86_dgab;
      od.df_dgamma_aa += 0.5*od.df_dgamma_ab;
      od.df_dgamma_bb += 0.5*od.df_dgamma_ab;
   }   
  
}

/////////////////////////////////////////////////////////////////////////////
// PBECFunctional

static ClassDesc PBECFunctional_cd(
  typeid(PBECFunctional),"PBECFunctional",1,"public DenFunctional",
  0, create<PBECFunctional>, create<PBECFunctional>);

PBECFunctional::PBECFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  local_ << SavableState::restore_state(s);
  s.get(gamma);
  s.get(beta);
}

PBECFunctional::PBECFunctional()
{
  local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);

  init_constants();
}

PBECFunctional::PBECFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  local_ << keyval->describedclassvalue("local");
  if (local_.null()) local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);

  init_constants();
  gamma = keyval->doublevalue("gamma", KeyValValuedouble(gamma));
  beta  = keyval->doublevalue("beta", KeyValValuedouble(beta));
}

void
PBECFunctional::init_constants()
{
  // in paper
  // gamma = 0.031091
  // beta  = 0.066725
  // in PBE.f:
  gamma = 0.03109069086965489503494086371273;
  beta  = 0.06672455060314922;
}

PBECFunctional::~PBECFunctional()
{
}

void
PBECFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  SavableState::save_state(local_.pointer(),s);
  s.put(gamma);
  s.put(beta);
}

int
PBECFunctional::need_density_gradient()
{
  return 1;
}

void
PBECFunctional::set_spin_polarized(int a)
{
  spin_polarized_ = a;
  local_->set_spin_polarized(a);
}

double
PBECFunctional::rho_deriv(double rho_a, double rho_b, double mdr,
                          double ec_local, double ec_local_dra)
{
  if (mdr < MIN_SQRTGAMMA) {
      return 0;
    }

  if (rho_b < MIN_DENSITY) {
#define Log(x) log(x)
#define Power(x,y) pow(x,y)
#define Pi M_PI
#define E M_E
      double ec = ec_local;
      double decdrhoa = ec_local_dra;
      double rhoa = rho_a;
      double result =
   (gamma*((-2*beta*Power(mdr,2)*Power(Pi,0.3333333333333333)*rhoa*
           (Power(beta,2)*decdrhoa*Power(E,(4*ec)/gamma)*Power(mdr,4)*
              Power(Pi,0.6666666666666666)*
              (Power(6,0.6666666666666666)*beta*Power(E,(2*ec)/gamma)*
                 Power(mdr,2)*Power(Pi,0.3333333333333333) - 
                96*(-1 + Power(E,(2*ec)/gamma))*gamma*
                 Power(rhoa,2.3333333333333335)) + 
             896*Power(6,0.3333333333333333)*
              Power(-1 + Power(E,(2*ec)/gamma),3)*Power(gamma,3)*
              Power(rhoa,3.6666666666666665)*
              (-(beta*Power(E,(2*ec)/gamma)*Power(mdr,2)*
                   Power(Pi,0.3333333333333333)) + 
                4*Power(6,0.3333333333333333)*(-1 + Power(E,(2*ec)/gamma))*
                 gamma*Power(rhoa,2.3333333333333335))))/
         (gamma*(Power(6,0.3333333333333333)*Power(beta,2)*
              Power(E,(2*ec)/gamma)*Power(mdr,4)*Power(Pi,0.6666666666666666)\
              - 8*Power(6,0.6666666666666666)*beta*
              (-1 + Power(E,(2*ec)/gamma))*gamma*Power(mdr,2)*
              Power(Pi,0.3333333333333333)*Power(rhoa,2.3333333333333335) + 
             384*Power(-1 + Power(E,(2*ec)/gamma),2)*Power(gamma,2)*
              Power(rhoa,4.666666666666667))*
           (Power(6,0.3333333333333333)*Power(beta,2)*Power(E,(4*ec)/gamma)*
              Power(mdr,4)*Power(Pi,0.6666666666666666) - 
             8*Power(6,0.6666666666666666)*beta*Power(E,(2*ec)/gamma)*
              (-1 + Power(E,(2*ec)/gamma))*gamma*Power(mdr,2)*
              Power(Pi,0.3333333333333333)*Power(rhoa,2.3333333333333335) + 
             384*Power(-1 + Power(E,(2*ec)/gamma),2)*Power(gamma,2)*
              Power(rhoa,4.666666666666667))) + 
        Log(1 + (beta*(-1 + Power(E,(2*ec)/gamma))*Power(mdr,2)*
             Power(Pi/6.,0.3333333333333333)*
             (Power(6,0.6666666666666666)*beta*Power(E,(2*ec)/gamma)*
                Power(mdr,2)*Power(Pi,0.3333333333333333) - 
               48*(-1 + Power(E,(2*ec)/gamma))*gamma*
                Power(rhoa,2.3333333333333335)))/
           (-(Power(6,0.3333333333333333)*Power(beta,2)*Power(E,(4*ec)/gamma)*
                Power(mdr,4)*Power(Pi,0.6666666666666666)) + 
             8*Power(6,0.6666666666666666)*beta*Power(E,(2*ec)/gamma)*
              (-1 + Power(E,(2*ec)/gamma))*gamma*Power(mdr,2)*
              Power(Pi,0.3333333333333333)*Power(rhoa,2.3333333333333335) - 
             384*Power(-1 + Power(E,(2*ec)/gamma),2)*Power(gamma,2)*
            Power(rhoa,4.666666666666667)))))/2.;
      return result;
    }
  if (rho_a < MIN_DENSITY) {
      // df_drho_a diverges for this case
      return 0.0;
    }

  double ra0,ra1,ra2,ra3,ra4,ra5,ra6,ra7,ra8,ra9,ra10,ra11,ra12,ra13,ra14,
      ra15,ra16,ra17,ra18,ra19,ra20,ra21,ra22,ra23,ra24,ra25,ra26,ra27,ra28,
      ra29,ra30,ra31,ra32,ra33,ra34,ra35;
  double dpbec_drho_a;

  ra0 = rho_b+rho_a;
  ra1 = -(1.0*rho_b)+rho_a;
  ra2 = 1/pow(ra0,2.0);
  ra3 = 1/pow(ra0,1.0);
  ra4 = -(1.0*ra1*ra3)+1.0;
  ra5 = ra1*ra3+1.0;
  ra6 = (0.66666666666666663*(ra3-(1.0*ra1*ra2)))/pow(ra5,
   0.33333333333333331)+(0.66666666666666663*(-(1.0*ra3)+ra1*ra2
   ))/pow(ra4,0.33333333333333331);
  ra7 = pow(ra5,0.66666666666666663)+pow(ra4,0.66666666666666663);
  ra8 = pow(ra7,2.0);
  ra9 = pow(mdr,2.0);
  ra10 = 1/pow(ra0,2.3333333333333335);
  ra11 = 1/ra8;
  ra12 = pow(ra7,3.0);
  ra13 = 1/ra12;
  ra14 = ec_local;
  ra15 = 1/pow(gamma,1.0);
  ra16 = 1/exp(8.0*ra13*ra14*ra15);
  ra17 = ra16-1.0;
  ra18 = 1/pow(ra17,1.0);
  ra19 = 0.25387282439081477*beta*ra9*ra10*ra11*ra18*ra15;
  ra20 = ra19+1.0;
  ra21 = pow(beta,2.0);
  ra22 = pow(mdr,4.0);
  ra23 = 1/pow(ra0,4.666666666666667);
  ra24 = 1/pow(ra7,4.0);
  ra25 = 1/pow(ra17,2.0);
  ra26 = 1/pow(gamma,2.0);
  ra27 = ra19+0.064451410964169495*ra21*ra22*ra23*ra24*ra25*ra26+
   1.0;
  ra28 = 1/pow(ra27,1.0);
  ra29 = 0.25387282439081477*beta*ra9*ra10*ra11*ra20*ra28*ra15+1.0
   ;
  ra30 = log(ra29);
  ra31 = 1/pow(ra0,3.3333333333333335);
  ra32 = -(0.50774564878162953*beta*ra9*ra10*ra6*ra13*ra18*ra15);
  ra33 = -(0.59236992357856788*beta*ra9*ra31*ra11*ra18*ra15);
  ra34 = -(8.0*ra13*ec_local_dra*ra15)+24.0*ra6*
   ra24*ra14*ra15;
  ra35 = -(0.25387282439081477*beta*ra9*ra10*ra11*ra25*ra16*ra34*
   ra15);
  dpbec_drho_a = (0.125*ra0*ra12*(-((0.25387282439081477*beta*ra9*
   ra10*ra11*ra20*(ra35+ra33+ra32-((0.12890282192833899*ra21*ra22*
   ra23*ra24*ra16*ra34*ra26)/pow(ra17,3.0))-((0.30077325116612436*
   ra21*ra22*ra24*ra25*ra26)/pow(ra0,5.666666666666667))-((
   0.25780564385667798*ra21*ra22*ra23*ra6*ra25*ra26)/pow(ra7,5.0))
   )*ra15)/pow(ra27,2.0))+0.25387282439081477*beta*ra9*ra10*ra11*
   ra28*(ra35+ra33+ra32)*ra15-(0.59236992357856788*beta*ra9*ra31*
   ra11*ra20*ra28*ra15)-(0.50774564878162953*beta*ra9*ra10*ra6*ra13*
   ra20*ra28*ra15))*gamma)/pow(ra29,1.0)+0.125*ra12*ra30*gamma+
   0.375*ra0*ra6*ra8*ra30*gamma;

  return dpbec_drho_a;
}

double
PBECFunctional::gab_deriv(double rho, double phi, double mdr, double ec_local)
{
  if (rho < MIN_DENSITY) return 0;

  if (mdr < MIN_SQRTGAMMA) {
      double result
          = (beta*phi*pow(M_PI/3.,1./3.))/(8.*pow(rho,4./3.));
      return result;
    }

  double g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16;
  double dpbec_dmdr;
  g0 = pow(phi,3.0);
  g1 = pow(beta,2.0);
  g2 = pow(mdr,3.0);
  g3 = 1/pow(phi,4.0);
  g4 = 1/pow(rho,4.666666666666667);
  g5 = 1/pow(gamma,1.0);
  g6 = 1/exp(1.0*1.0/g0*ec_local*g5)-1.0;
  g7 = 1/pow(g6,1.0);
  g8 = 1/pow(g6,2.0);
  g9 = 1/pow(gamma,2.0);
  g10 = pow(mdr,2.0);
  g11 = 1/pow(phi,2.0);
  g12 = 1/pow(rho,2.3333333333333335);
  g13 = 0.063468206097703692*beta*g10*g11*g12*g7*g5;
  g14 = g13+0.0040282131852605934*g1*pow(mdr,4.0)*g3*g4*g8*g9+
   1.0;
  g15 = 1/pow(g14,1.0);
  g16 = g13+1.0;
  dpbec_dmdr = (g0*rho*(0.12693641219540738*beta*mdr*g11*g12*g16*g15
   *g5-((0.063468206097703692*beta*g10*g11*g12*(
   0.12693641219540738*beta*mdr*g11*g12*g7*g5+0.016112852741042374
   *g1*g2*g3*g4*g8*g9)*g16*g5)/pow(g14,2.0))+0.0080564263705211869
   *g1*g2*g3*g4*g7*g15*g9)*gamma)/pow(0.063468206097703692*beta*g10*
   g11*g12*g16*g15*g5+1.0,1.0);

  return dpbec_dmdr/mdr;
}

void
PBECFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  double ec_local, dec_local_rs, dec_local_zeta;
  local_->point_lc(id, od, ec_local, dec_local_rs, dec_local_zeta); 
  double rho = id.a.rho+id.b.rho;
  double rho_13 = pow(rho, 1./3.);
  double rho_43 = rho*rho_13;
  double zeta = (id.a.rho - id.b.rho)/rho;
  double phi = 0.5*(pow(1+zeta, 2./3.)+pow(1-zeta, 2./3.));
  double mdr = sqrt(id.a.gamma + id.b.gamma + 2*id.gamma_ab);

  double pbec;

  double e0,e1,e2,e3,e4,e5,e6;
  e0 = pow(phi,3.0);
  e1 = pow(mdr,2.0);
  e2 = 1/pow(phi,2.0);
  e3 = 1/pow(rho,2.3333333333333335);
  e4 = 1/pow(gamma,1.0);
  e5 = 1/exp(1.0*1.0/e0*ec_local*e4)-1.0;
  e6 = (0.063468206097703692*beta*e1*e2*e3*e4)/pow(e5,1.0);
  pbec = e0*rho*log((0.063468206097703692*beta*e1*e2*e3*(e6+1.0)*
     e4)/pow(e6+(0.0040282131852605934*pow(beta,2.0)*pow(mdr,4.0))
     /pow(phi,4.0)/pow(rho,4.666666666666667)/pow(e5,2.0)/pow(
     gamma,2.0)+1.0,1.0)+1.0)*gamma;

  if (compute_potential_) {
      double drs_drho_a = -0.20678349696646667/rho_43; // == drs_drho_b
      double dzeta_drho_a = 0.;
      if (zeta < MAX_ZETA) dzeta_drho_a = 1./rho * ( 1. - zeta);
      double dzeta_drho_b = 0.;
      if (zeta > MIN_ZETA) dzeta_drho_b = 1./rho * (-1. - zeta);

      //double ec_local_dra = od.df_drho_a;
      //double ec_local_drb = od.df_drho_b;
      double ec_local_dra
          = dec_local_rs*drs_drho_a + dec_local_zeta*dzeta_drho_a;
      double ec_local_drb
          = dec_local_rs*drs_drho_a + dec_local_zeta*dzeta_drho_b;

      double df_drho_a = rho_deriv(id.a.rho, id.b.rho, mdr,
                                   ec_local, ec_local_dra);
      double df_drho_b = rho_deriv(id.b.rho, id.a.rho, mdr,
                                   ec_local, ec_local_drb);

      double df_dgab = gab_deriv(rho, phi, mdr, ec_local);

      od.df_drho_a += df_drho_a;
      od.df_drho_b += df_drho_b;
      od.df_dgamma_aa += 0.5*df_dgab;
      od.df_dgamma_bb += 0.5*df_dgab;
      od.df_dgamma_ab += df_dgab;
    }

  od.energy += pbec;

//   cout << scprintf("id.a.del_rho = %12.8f %12.8f %12.8f", id.a.del_rho[0],
//                    id.a.del_rho[1], id.a.del_rho[2])
//        << std::endl;
//   cout << scprintf("id.b.del_rho = %12.8f %12.8f %12.8f", id.b.del_rho[0],
//                    id.b.del_rho[1], id.b.del_rho[2])
//        << std::endl;
//   cout << "id.a.gamma = " << id.a.gamma << std::endl;
//   cout << "id.b.gamma = " << id.b.gamma << std::endl;
//   cout << "od.df_drho_a = " << od.df_drho_a << std::endl;
//   cout << "od.df_drho_b = " << od.df_drho_b << std::endl;
//   cout << "od.df_dgamma_aa = " << od.df_dgamma_aa << std::endl;
//   cout << "od.df_dgamma_ab = " << od.df_dgamma_ab << std::endl;
//   cout << "od.df_dgamma_bb = " << od.df_dgamma_bb << std::endl;
}

/////////////////////////////////////////////////////////////////////////////
// PW91CFunctional

static ClassDesc PW91CFunctional_cd(
  typeid(PW91CFunctional),"PW91CFunctional",1,"public DenFunctional",
  0, create<PW91CFunctional>, create<PW91CFunctional>);

PW91CFunctional::PW91CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  local_ << SavableState::restore_state(s);
  init_constants();
}

PW91CFunctional::PW91CFunctional()
{
  local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);
  init_constants();
}

PW91CFunctional::PW91CFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  local_ << keyval->describedclassvalue("local");
  if (local_.null()) local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);
  init_constants();
}

void
PW91CFunctional::init_constants()
{
  a = 23.266;
  b = 7.389e-3;
  c = 8.723;
  d = 0.472;
  alpha = 0.09;
  c_c0 = 0.004235;
  // c_x = -0.001667 in  Phys Rev B 46, p 6671, Perdew, et. al.
  // c_x as in PBE.f:
  c_x = -0.001667212;
}

PW91CFunctional::~PW91CFunctional()
{
}

void
PW91CFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  SavableState::save_state(local_.pointer(),s);
}

int
PW91CFunctional::need_density_gradient()
{
  return 1;
}

void
PW91CFunctional::set_spin_polarized(int a)
{
  spin_polarized_ = a;
  local_->set_spin_polarized(a);
}

double
PW91CFunctional::limit_df_drhoa(double rhoa, double mdr,
                                double ec, double decdrhoa)
{
  double result;
  // blame mathematica
  double v = 15.755920349483144;
  double cx = c_x;
  double cc0 = c_c0;
  double beta = v * cc0;
  double e2ec = Power(E,(2*alpha*ec)/Power(beta,2));
  double e4ec = e2ec * e2ec;
  double e8ec = e4ec * e4ec;
  double e25mdr2 = Power(E,-(25*Power(mdr,2))/
         (Power(6,0.6666666666666666)*Power(Pi,1.3333333333333333)*
          Power(rhoa,2.6666666666666665)));
  result =
   (Power(mdr,2)*Power(Pi/6.,0.3333333333333333)*
       (1750*Power(6,0.3333333333333333)*a*Power(Pi,0.6666666666666666) - 
         1750000*Power(6,0.3333333333333333)*c*cc0*
          Power(Pi,0.6666666666666666) - 
         2500000*Power(6,0.3333333333333333)*c*cx*
          Power(Pi,0.6666666666666666) + 
         (8988*Pi)/Power(1/rhoa,0.3333333333333333) - 
         (3500000*cc0*Pi)/Power(1/rhoa,0.3333333333333333) - 
         (5000000*cx*Pi)/Power(1/rhoa,0.3333333333333333) + 
         875*Power(6,0.6666666666666666)*b*Power(Pi,0.3333333333333333)*
          Power(1/rhoa,0.3333333333333333) - 
         875000*Power(6,0.6666666666666666)*cc0*d*
          Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
         1250000*Power(6,0.6666666666666666)*cx*d*
          Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
         26250000*b*cc0*Power(1/rhoa,0.6666666666666666) - 
         37500000*b*cx*Power(1/rhoa,0.6666666666666666))*v)
     *e25mdr2/(1.4e7*
       (2*Power(6,0.3333333333333333)*c*Power(Pi,0.6666666666666666) + 
         (4*Pi)/Power(1/rhoa,0.3333333333333333) + 
         Power(6,0.6666666666666666)*d*Power(Pi,0.3333333333333333)*
          Power(1/rhoa,0.3333333333333333) + 
         30*b*Power(1/rhoa,0.6666666666666666))*Power(rhoa,2.3333333333333335)
         );
  result +=
          rhoa*((Power(beta,2)*
          ((-2*alpha*(-(alpha*e4ec*Power(mdr,4)*
                     Power(Pi/6.,0.6666666666666666))/
                  (32.*beta*(-1 + e4ec)*Power(rhoa,4.666666666666667)) + 
                 (Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                  (8.*Power(rhoa,2.3333333333333335)))*
               ((-7*Power(alpha,2)*e8ec*Power(mdr,4)*
                    Power(Pi/6.,0.6666666666666666))/
                  (24.*beta*Power(-1 + e4ec,2)*Power(rhoa,5.666666666666667))\
                  - (Power(alpha,3)*decdrhoa*e8ec*Power(mdr,4)*
                    Power(Pi/6.,0.6666666666666666))/
                  (2.*Power(beta,3)*Power(-1 + e4ec,3)*
                    Power(rhoa,4.666666666666667)) + 
                 (7*alpha*e4ec*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                  (12.*(-1 + e4ec)*Power(rhoa,3.3333333333333335)) + 
                 (Power(alpha,2)*decdrhoa*e4ec*Power(mdr,2)*
                    Power(Pi/6.,0.3333333333333333))/
                  (Power(beta,2)*Power(-1 + e4ec,2)*
                    Power(rhoa,2.3333333333333335))))/
             Power(beta + (Power(alpha,2)*e8ec*Power(mdr,4)*
                  Power(Pi/6.,0.6666666666666666))/
                (16.*beta*Power(-1 + e4ec,2)*Power(rhoa,4.666666666666667)) - 
               (alpha*e4ec*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                (4.*(-1 + e4ec)*Power(rhoa,2.3333333333333335)),2) + 
            (2*alpha*((7*alpha*e4ec*Power(mdr,4)*
                    Power(Pi/6.,0.6666666666666666))/
                  (48.*beta*(-1 + e4ec)*Power(rhoa,5.666666666666667)) + 
                 (Power(alpha,2)*decdrhoa*e4ec*Power(mdr,4)*
                    Power(Pi/6.,0.6666666666666666))/
                  (8.*Power(beta,3)*Power(-1 + e4ec,2)*
                    Power(rhoa,4.666666666666667)) - 
                 (7*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                  (24.*Power(rhoa,3.3333333333333335))))/
             (beta + (Power(alpha,2)*e8ec*Power(mdr,4)*
                  Power(Pi/6.,0.6666666666666666))/
                (16.*beta*Power(-1 + e4ec,2)*Power(rhoa,4.666666666666667)) - 
               (alpha*e4ec*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                (4.*(-1 + e4ec)*Power(rhoa,2.3333333333333335)))))/
        (4.*alpha*(1 + (2*alpha*
               (-(alpha*e4ec*Power(mdr,4)*Power(Pi/6.,0.6666666666666666))/
                  (32.*beta*(-1 + e4ec)*Power(rhoa,4.666666666666667)) + 
                 (Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                  (8.*Power(rhoa,2.3333333333333335))))/
             (beta + (Power(alpha,2)*e8ec*Power(mdr,4)*
                  Power(Pi/6.,0.6666666666666666))/
                (16.*beta*Power(-1 + e4ec,2)*Power(rhoa,4.666666666666667)) - 
               (alpha*e4ec*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
                (4.*(-1 + e4ec)*Power(rhoa,2.3333333333333335))))) + 
       (Power(mdr,4)*(1750*Power(6,0.3333333333333333)*a*
             Power(Pi,0.6666666666666666) - 
            1750000*Power(6,0.3333333333333333)*c*cc0*
             Power(Pi,0.6666666666666666) - 
            2500000*Power(6,0.3333333333333333)*c*cx*
             Power(Pi,0.6666666666666666) + 
            (8988*Pi)/Power(1/rhoa,0.3333333333333333) - 
            (3500000*cc0*Pi)/Power(1/rhoa,0.3333333333333333) - 
            (5000000*cx*Pi)/Power(1/rhoa,0.3333333333333333) + 
            875*Power(6,0.6666666666666666)*b*Power(Pi,0.3333333333333333)*
             Power(1/rhoa,0.3333333333333333) - 
            875000*Power(6,0.6666666666666666)*cc0*d*
             Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
            1250000*Power(6,0.6666666666666666)*cx*d*
             Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
            26250000*b*cc0*Power(1/rhoa,0.6666666666666666) - 
            37500000*b*cx*Power(1/rhoa,0.6666666666666666))*v)
        *e25mdr2/(1.26e6*Pi*
          (2*Power(6,0.3333333333333333)*c*Power(Pi,0.6666666666666666) + 
            (4*Pi)/Power(1/rhoa,0.3333333333333333) + 
            Power(6,0.6666666666666666)*d*Power(Pi,0.3333333333333333)*
             Power(1/rhoa,0.3333333333333333) + 
            30*b*Power(1/rhoa,0.6666666666666666))*Power(rhoa,6)) - 
       (Power(mdr,2)*Power(Pi/6.,0.3333333333333333)*
          (1750*Power(6,0.3333333333333333)*a*Power(Pi,0.6666666666666666) - 
            1750000*Power(6,0.3333333333333333)*c*cc0*
             Power(Pi,0.6666666666666666) - 
            2500000*Power(6,0.3333333333333333)*c*cx*
             Power(Pi,0.6666666666666666) + 
            (8988*Pi)/Power(1/rhoa,0.3333333333333333) - 
            (3500000*cc0*Pi)/Power(1/rhoa,0.3333333333333333) - 
            (5000000*cx*Pi)/Power(1/rhoa,0.3333333333333333) + 
            875*Power(6,0.6666666666666666)*b*Power(Pi,0.3333333333333333)*
             Power(1/rhoa,0.3333333333333333) - 
            875000*Power(6,0.6666666666666666)*cc0*d*
             Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
            1250000*Power(6,0.6666666666666666)*cx*d*
             Power(Pi,0.3333333333333333)*Power(1/rhoa,0.3333333333333333) - 
            26250000*b*cc0*Power(1/rhoa,0.6666666666666666) - 
            37500000*b*cx*Power(1/rhoa,0.6666666666666666))*v)*e25mdr2/
        (6.e6*
          (2*Power(6,0.3333333333333333)*c*Power(Pi,0.6666666666666666) + 
            (4*Pi)/Power(1/rhoa,0.3333333333333333) + 
            Power(6,0.6666666666666666)*d*Power(Pi,0.3333333333333333)*
             Power(1/rhoa,0.3333333333333333) + 
            30*b*Power(1/rhoa,0.6666666666666666))*
          Power(rhoa,3.3333333333333335)) + 
       (Power(mdr,2)*Power(Pi/6.,0.3333333333333333)*
          (1875*Power(6,0.6666666666666666)*Power(b,2)*
             Power(Pi,0.3333333333333333) - 
            (500*Power(6,0.3333333333333333)*a*Power(Pi,1.6666666666666667))/
             Power(1/rhoa,1.3333333333333333) + 
            (1284*Power(6,0.3333333333333333)*c*Power(Pi,1.6666666666666667))/
             Power(1/rhoa,1.3333333333333333) + 
            (57780*b*Pi)/Power(1/rhoa,0.6666666666666666) - 
            (750*b*c*Pi)/Power(1/rhoa,0.6666666666666666) + 
            (750*a*d*Pi)/Power(1/rhoa,0.6666666666666666) + 
            (7500*Power(6,0.3333333333333333)*a*b*
               Power(Pi,0.6666666666666666))/Power(1/rhoa,0.3333333333333333)\
             - 500*Power(6,0.6666666666666666)*b*Power(Pi,1.3333333333333333)*
             rhoa + 1284*Power(6,0.6666666666666666)*d*
             Power(Pi,1.3333333333333333)*rhoa)*v)*e25mdr2/
        (3.e6*
          Power(2*Power(6,0.3333333333333333)*c*
             Power(Pi,0.6666666666666666) + 
            (4*Pi)/Power(1/rhoa,0.3333333333333333) + 
            Power(6,0.6666666666666666)*d*Power(Pi,0.3333333333333333)*
             Power(1/rhoa,0.3333333333333333) + 
            30*b*Power(1/rhoa,0.6666666666666666),2)*
          Power(rhoa,4.333333333333333))) + 
    (Power(beta,2)*Log(1 + (2*alpha*
            (-(alpha*e4ec*Power(mdr,4)*Power(Pi/6.,0.6666666666666666))/
               (32.*beta*(-1 + e4ec)*Power(rhoa,4.666666666666667)) + 
              (Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
               (8.*Power(rhoa,2.3333333333333335))))/
          (beta + (Power(alpha,2)*e8ec*Power(mdr,4)*
               Power(Pi/6.,0.6666666666666666))/
             (16.*beta*Power(-1 + e4ec,2)*Power(rhoa,4.666666666666667)) - 
            (alpha*e4ec*Power(mdr,2)*Power(Pi/6.,0.3333333333333333))/
             (4.*(-1 + e4ec)*Power(rhoa,2.3333333333333335)))))/(4.*alpha)
          ;

  return result;
}

void
PW91CFunctional::point(const PointInputData &id, PointOutputData &od)
{
  double ec_local, dec_local_rs, dec_local_zeta;
  local_->point_lc(id, od, ec_local, dec_local_rs, dec_local_zeta); 
  double rho = id.a.rho+id.b.rho;
  double rho_13 = pow(rho, 1./3.);
  double rho_43 = rho*rho_13;
  double rs = 0.62035049089940009/rho_13;
  double z = (id.a.rho - id.b.rho)/rho;
  double gamma = sqrt(id.a.gamma + id.b.gamma + 2*id.gamma_ab);

  double pwc, dpwc_drs, dpwc_dg;

  if (rho < MIN_DENSITY) return;
  if (gamma < MIN_SQRTGAMMA) {
      if (!compute_potential_) return;

      double limit_dpwc_dgaa, limit_dpwc_dgab;

      double ga0;
      ga0 = 6.2337093539550051e7*b*c_x*pow(rs,7.0)+(6233709.3539550053
       *c_x*d-(4363.596547768504*b))*pow(rs,6.0)+(6233709.3539550053
       *c*c_x-(4363.596547768504*a))*pow(rs,5.0)+(6233709.3539550053
       *c_x-11205.715934669519)*pow(rs,4.0);
      limit_dpwc_dgaa = -((1.0*(1.4645918875615231*ga0*pow(z+1.0,
       0.66666666666666663)+1.4645918875615231*ga0*pow(-(1.0*z)+
       1.0,0.66666666666666663)))/pow(1.8929525610284735e7*b*pow(rs,
       3.0)+1892952.5610284733*d*pow(rs,2.0)+1892952.5610284733*c*
       rs+1892952.5610284733,1.0));

      double gab0;
      gab0 = 6.2337093539550051e7*b*c_x*pow(rs,7.0)+(
       6233709.3539550053*c_x*d-(4363.596547768504*b))*pow(rs,6.0)+(
       6233709.3539550053*c*c_x-(4363.596547768504*a))*pow(rs,5.0)+(
       6233709.3539550053*c_x-11205.715934669519)*pow(rs,4.0);
      limit_dpwc_dgab = -((1.0*(1.4645918875615231*gab0*pow(z+1.0,
       0.66666666666666663)+1.4645918875615231*gab0*pow(-(1.0*z)+
       1.0,0.66666666666666663)))/pow(9464762.8051423673*b*pow(rs,
       3.0)+946476.28051423666*d*pow(rs,2.0)+946476.28051423666*c*
       rs+946476.28051423666,1.0));

      od.df_dgamma_aa += limit_dpwc_dgaa;
      od.df_dgamma_bb += limit_dpwc_dgaa;
      od.df_dgamma_ab += limit_dpwc_dgab;

      return;
    }

  double e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15;
  e2 = rs*rs; // pow(rs,2.0);
  e0 = e2*rs; // pow(rs,3.0);
  e1 = e0*e0*rs; // pow(rs,7.0);
  e3 = pow(z+1.0,2./3.)+pow(-z+1.0, 2./3.);
  e4 = gamma*gamma; // pow(gamma,2.0);
  e5 = e3*e3; // pow(e3,2.0);
  e6 = c_c0*c_c0; // pow(c_c0,2.0);
  e7 = e3*e3*e3; // pow(e3,3.0);
  e8 = 1./c_c0; // 1/pow(c_c0,1.0);
  e9 = 1/e5;
  e10 = 1/e6;
  e11 = exp(-0.064451410964169495*alpha*e10*ec_local/e7)-1.0;
  e12 = 1./e11; // 1/pow(e11,1.0);
  e13 = e1*e1; // pow(rs,14.0);
  e14 = 1/(e7*e3); // pow(e3,4.0);
  e15 = e4*e4; // pow(gamma,4.0);
  pwc = (0.238732414637843*((15.515564128703568*e6*e7*log((
     0.12693641219540738*alpha*e8*(6.5448368524534253*alpha*e8*e13*
     e14*e12*e15+7.18052672676657*e1*e9*e4))/((
     0.83077810845472067*alpha*alpha*e10*e13*e14*e15)/(e11*e11)
     +0.91147030036898102*alpha*e8*e1*e9*e12*e4+1.0)+
     1.0))/alpha+(14.141975896783627*e1*((0.001*(b*e2+a
     *rs+2.568))/(10.0*b*e0+d*e2+c*rs+1.0)-(
     1.4285714285714286*c_x)-(1.0*c_c0))*e3*e4)*exp(-
     29.773894288156065*e1*rs*e5*e4)))/e0;

  if (compute_potential_) {

      double drs_drhoa = -0.20678349696646667/rho_43;
      double drs_drhob = drs_drhoa;

      double r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18
          ,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,t0;

      r0 = pow(rs,3.0);
      r1 = 1/pow(alpha,1.0);
      r2 = pow(c_c0,2.0);
      r3 = pow(z+1.0,0.66666666666666663)+pow(-(1.0*z)+1.0,
                                                  0.66666666666666663);
      r4 = pow(r3,3.0);
      r5 = 1/pow(c_c0,1.0);
      r6 = pow(rs,7.0);
      r7 = pow(r3,2.0);
      r8 = 1/r7;
      r9 = 1/r2;
      r10 = exp(-0.064451410964169495*alpha*r9*ec_local*1.0/r4);
      r11 = r10-1.0;
      r12 = 1/pow(r11,1.0);
      r13 = pow(gamma,2.0);
      r14 = pow(alpha,2.0);
      r15 = pow(rs,14.0);
      r16 = 1/pow(r3,4.0);
      r17 = 1/pow(r11,2.0);
      r18 = pow(gamma,4.0);
      r19 = 0.83077810845472067*r14*r9*r15*r16*r17*r18+
            0.91147030036898102*alpha*r5*r6*r8*r12*r13+1.0;
      r20 = 1/pow(r19,1.0);
      r21 = 6.5448368524534253*alpha*r5*r15*r16*r12*r18+
            7.18052672676657*r6*r8*r13;
      r22 = 0.12693641219540738*alpha*r5*r20*r21+1.0;
      r23 = pow(rs,6.0);
      r24 = 1/pow(c_c0,3.0);
      r25 = dec_local_rs;
      r26 = pow(rs,13.0);
      r27 = 1/pow(r3,7.0);
      r28 = pow(rs,2.0);
      r29 = b*r28+a*rs+2.568;
      r30 = 10.0*b*r0+d*r28+c*rs+1.0;
      r31 = 1/pow(r30,1.0);
      r32 = exp(-29.773894288156065*pow(rs,8.0)*r7*r13);
      r33 = 0.001*r29*r31-(1.4285714285714286*c_x)-(1.0*c_c0);
      t0 = -((0.71619724391352901*(15.515564128703568*r1*r2*r4*log
         (r22)+14.141975896783627*r6*r33*r3*r13*r32))/pow(rs,4.0));
      dpwc_drs = t0+(0.238732414637843*(-(3368.493563011893*r15*
         r33*r4*r18*r32)+98.993831277485398*r23*r33*r3*r13*r32+
         14.141975896783627*r6*(0.001*(2.0*b*rs+a)*r31-((0.001*
         r29*(30.0*b*r28+2.0*d*rs+c))/pow(r30,2.0)))*r3*r13*r32+(
         15.515564128703568*r1*r2*r4*(0.12693641219540738*alpha*r5*
         r20*(0.42182396967091729*r14*r24*r25*r15*r27*r17*r10*r18+
         91.627715934347961*alpha*r5*r26*r16*r12*r18+
         50.263687087365994*r23*r8*r13)-((0.12693641219540738*alpha*
         r5*r21*((0.10708964257610123*pow(alpha,3.0)*r25*r15*r27*r10
         *r18)/pow(c_c0,4.0)/pow(r11,3.0)+11.630893518366092*r14*
         r9*r26*r16*r17*r18+(0.058745546910716179*r14*r24*r25*r6*r17*
         r10*r13)/pow(r3,5.0)+6.3802921025828665*alpha*r5*r23*r8*r12
         *r13))/pow(r19,2.0))))/pow(r22,1.0)))/r0;

      if (id.a.rho > MIN_DENSITY && id.b.rho > MIN_DENSITY) {
          double dpwc_dz;

          double z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17
              ,z18 ,z19,z20,z21,z22,z23,z24,z25,z26,z27,z28,z29,z30,z31;
          z0 = pow(rs,3.0);
          z1 = 1/pow(alpha,1.0);
          z2 = pow(c_c0,2.0);
          z3 = -(1.0*z)+1.0;
          z4 = z+1.0;
          z5 = pow(z4,0.66666666666666663)+pow(z3,0.66666666666666663);
          z6 = pow(z5,3.0);
          z7 = 1/pow(c_c0,1.0);
          z8 = pow(rs,7.0);
          z9 = pow(z5,2.0);
          z10 = 1/z9;
          z11 = 1/z2;
          z12 = 1/z6;
          z13 = exp(-0.064451410964169495*alpha*z11*ec_local*z12);
          z14 = z13-1.0;
          z15 = 1/pow(z14,1.0);
          z16 = pow(gamma,2.0);
          z17 = pow(alpha,2.0);
          z18 = pow(rs,14.0);
          z19 = 1/pow(z5,4.0);
          z20 = 1/pow(z14,2.0);
          z21 = pow(gamma,4.0);
          z22 = 0.83077810845472067*z17*z11*z18*z19*z20*z21+
                0.91147030036898102*alpha*z7*z8*z10*z15*z16+1.0;
          z23 = 1/pow(z22,1.0);
          z24 = 6.5448368524534253*alpha*z7*z18*z19*z15*z21+
                7.18052672676657*z8*z10*z16;
          z25 = 0.12693641219540738*alpha*z7*z23*z24+1.0;
          z26 = 0.66666666666666663/pow(z4,0.33333333333333331)-(
              0.66666666666666663/pow(z3,0.33333333333333331));
          z27 = -(0.064451410964169495*alpha*z11*dec_local_zeta*z12)
               +0.19335423289250847*alpha*z11*ec_local*z26*z19;
          z28 = 1/pow(z5,5.0);
          z29 = pow(rs,2.0);
          z30 = (0.001*(b*z29+a*rs+2.5680000000000001))/pow(10.0*b*z0+d*
              z29+c*rs+1.0,1.0)-(1.4285714285714286*c_x)-(1.0*c_c0);
          z31 = exp(-29.773894288156065*pow(rs,8.0)*z9*z16);
          dpwc_dz = (0.238732414637843*(46.546692386110699*z1*z2*z26*z9*
           log(z25)-(842.12339075297325*pow(rs,15.0)*z30*z26*z9*z21*z31)+
           14.141975896783627*z8*z30*z26*z16*z31+(15.515564128703568*z1*z2
           *z6*(0.12693641219540738*alpha*z7*z23*(-(6.5448368524534253*
           alpha*z7*z18*z19*z27*z20*z13*z21)-(26.179347409813701*alpha*z7*
           z18*z26*z28*z15*z21)-(14.36105345353314*z8*z26*z12*z16))-((
           0.12693641219540738*alpha*z7*z24*(-((1.6615562169094413*z17*z11
           *z18*z19*z27*z13*z21)/pow(z14,3.0))-(3.3231124338188827*z17*z11
           *z18*z26*z28*z20*z21)-(0.91147030036898102*alpha*z7*z8*z10*z27*
           z20*z13*z16)-(1.822940600737962*alpha*z7*z8*z26*z12*z15*z16)))/
           pow(z22,2.0))))/pow(z25,1.0)))/z0;

          double dz_drhoa = ( 1. - (id.a.rho-id.b.rho)/rho)/rho;
          double dz_drhob = (-1. - (id.a.rho-id.b.rho)/rho)/rho;

          od.df_drho_a += dpwc_drs * drs_drhoa + dpwc_dz * dz_drhoa;
          od.df_drho_b += dpwc_drs * drs_drhob + dpwc_dz * dz_drhob;
        }
      else if (id.a.rho > MIN_DENSITY) {
          // df_drho_b diverges
          double dzeta_drhoa = 1./rho * (1. - z);
          double drs_drhoa = -rs/(3.*rho);
          double dec_local_drhoa = dec_local_rs*drs_drhoa
                                 + dec_local_zeta*dzeta_drhoa;
          od.df_drho_a += limit_df_drhoa(id.a.rho, gamma,
                                         ec_local, dec_local_drhoa);
        }
      else if (id.b.rho > MIN_DENSITY) {
          // df_drho_a diverges
          double dzeta_drhob = 1./rho * (-1. - z);
          double drs_drhob = -rs/(3.*rho);
          double dec_local_drhob = dec_local_rs*drs_drhob
                                 + dec_local_zeta*dzeta_drhob;
          od.df_drho_b += limit_df_drhoa(id.b.rho, gamma,
                                         ec_local, dec_local_drhob);
        }

      double g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18
          ,g19,g20,g21,g22,g23;
      g0 = pow(rs,3.0);
      g1 = pow(c_c0,2.0);
      g2 = pow(z+1.0,0.66666666666666663)+pow(-(1.0*z)+1.0,
                                              0.66666666666666663);
      g3 = pow(g2,3.0);
      g4 = 1/pow(c_c0,1.0);
      g5 = pow(rs,7.0);
      g6 = pow(g2,2.0);
      g7 = 1/g6;
      g8 = pow(rs,14.0);
      g9 = 1/pow(g2,4.0);
      g10 = 1/g1;
      g11 = exp(-0.064451410964169495*alpha*g10*ec_local*1.0/g3)-
            1.0;
      g12 = 1/pow(g11,1.0);
      g13 = pow(gamma,3.0);
      g14 = pow(gamma,2.0);
      g15 = pow(alpha,2.0);
      g16 = 1/pow(g11,2.0);
      g17 = pow(gamma,4.0);
      g18 = 0.83077810845472067*g15*g10*g8*g9*g16*g17+
            0.91147030036898102*alpha*g4*g5*g7*g12*g14+1.0;
      g19 = 1/pow(g18,1.0);
      g20 = 6.5448368524534253*alpha*g4*g8*g9*g12*g17+7.18052672676657
           *g5*g7*g14;
      g21 = pow(rs,2.0);
      g22 = (0.001*(b*g21+a*rs+2.5680000000000001))/pow(10.0*b*g0+d*
            g21+c*rs+1.0,1.0)-(1.4285714285714286*c_x)-(1.0*c_c0);
      g23 = exp(-29.773894288156065*pow(rs,8.0)*g6*g14);
      dpwc_dg = (0.238732414637843*(-(842.12339075297325*pow(rs,15.0
       )*g22*g3*g13*g23)+28.283951793567255*g5*g22*g2*gamma*g23+(
       15.515564128703568*g1*g3*(-((0.12693641219540738*alpha*g4*(
       3.3231124338188827*g15*g10*g8*g9*g16*g13+1.822940600737962*
       alpha*g4*g5*g7*g12*gamma)*g20)/pow(g18,2.0))+
       0.12693641219540738*alpha*g4*(26.179347409813701*alpha*g4*g8*g9
       *g12*g13+14.36105345353314*g5*g7*gamma)*g19))/pow(alpha,1.0)/
       pow(0.12693641219540738*alpha*g4*g19*g20+1.0,1.0)))/g0;

      od.df_dgamma_aa += dpwc_dg / (2 * gamma);
      od.df_dgamma_bb += dpwc_dg / (2 * gamma);
      od.df_dgamma_ab += dpwc_dg / gamma;
    }

  od.energy += pwc;
}

/////////////////////////////////////////////////////////////////////////////
// PW91XFunctional

static ClassDesc PW91XFunctional_cd(
  typeid(PW91XFunctional),"PW91XFunctional",1,"public DenFunctional",
  0, create<PW91XFunctional>, create<PW91XFunctional>);

PW91XFunctional::PW91XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(a);
  s.get(b);
  s.get(c);
  s.get(d);
  s.get(a_x);
}

PW91XFunctional::PW91XFunctional()
{
  init_constants();
}

PW91XFunctional::PW91XFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  a_x = keyval->doublevalue("a_x", KeyValValuedouble(a_x));
  a = keyval->doublevalue("a", KeyValValuedouble(a));
  b = keyval->doublevalue("b", KeyValValuedouble(b));
  c = keyval->doublevalue("c", KeyValValuedouble(c));
  d = keyval->doublevalue("d", KeyValValuedouble(d));
}

PW91XFunctional::~PW91XFunctional()
{
}

void
PW91XFunctional::init_constants()
{
  a = 0.19645;
  b = 7.7956;
  c = 0.2743;
  // the PW91 paper had d = 0.1508
  // PBE.f has the following
  d = 0.15084;
  // a_x -(3/4)*(3/pi)^(1/3)
  // the following rounded a_x appears in PBE.f
  a_x = -0.7385588;
}

void
PW91XFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(a);
  s.put(b);
  s.put(c);
  s.put(d);
  s.put(a_x);
}

int
PW91XFunctional::need_density_gradient()
{
  return 1;
}

void
PW91XFunctional::spin_contrib(const PointInputData::SpinData &i,
                               double &pw, double &dpw_dr, double &dpw_dg)
{
  double rho = i.rho;
  double gamma = i.gamma;
  const double pi = M_PI;

  if (rho < MIN_DENSITY) {
      pw = 0.;
      dpw_dr = 0.;
      dpw_dg = 0.; // really -inf
      return;
    }
  if (gamma < MIN_GAMMA) {
      double rho_43 = rho * i.rho_13; // rho^(4/3)
      // 2^(1/3) * a_x * rho^(4/3)
      pw = 1.2599210498948732 * a_x * rho_43;
      // (4/3) 2^(1/3) a_x rho^(1/3)
      dpw_dr = 1.6798947331931642 * a_x * i.rho_13;
      // (6^(2/3) a_x / (24 3^(1/3) pi^4/3)) (d - c) / rho^(4/3)
      dpw_dg = - 0.020732388737701564 * a_x * (d - c) / rho_43;
      return;
    }

  // this has been generated by macsyma and modified
  double t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
  t0 = 1.8171205928321397; // pow(6,1.0/3.0);
  t1 = rho * i.rho_13; // pow(rho,4.0/3.0);
  t2 = 1/t0;
  t3 = 0.46619407703541166; // 1/pow(pi,2.0/3.0);
  t4 = 1/t1;
  t5 = sqrt(gamma);
  t6 = (t2*a*t3*t4*asinh((t2*b*t3*t4*t5)/2.0)*t5)/2.0;
  t7 = 0.30285343213868998; // 1/pow(6,2.0/3.0);
  t8 = 0.21733691746289932; // 1/pow(pi,4.0/3.0);
  t9 = t4*t4; // 1/pow(rho,8.0/3.0);
  pw = (t0*a_x*t1*((t7*t8*t9*gamma*(-(d*exp(-25.0*t7*t8*t9*gamma))+c)
     )/4.0+t6+1))/pow(3,1.0/3.0)/((t2*pow(gamma,2.0))/24000.0/pow(pi,8.0
     /3.0)/pow(rho,16.0/3.0)+t6+1);
  if (compute_potential_) {
      double r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18
          ,r19,r20,r21,r22;
      r0 = 0.69336127435063466; // 1/pow(3,1.0/3.0);
      r1 = t0; // pow(6,1.0/3.0);
      r2 = t1; // pow(rho,4.0/3.0);
      r3 = 1/r1;
      r4 = t3; // 1/pow(pi,2.0/3.0);
      r5 = 1/r2;
      r6 = t5; // sqrt(gamma);
      r7 = asinh((r3*b*r4*r5*r6)/2.0);
      r8 = (r3*a*r4*r5*r7*r6)/2.0;
      r9 = 1/pow(pi,8.0/3.0);
      r10 = gamma*gamma; // pow(gamma,2.0);
      r11 = (r3*r9*r10)/24000.0/pow(rho,16.0/3.0)+r8+1;
      r12 = 1/r11;
      r13 = -((2.0/3.0*r3*a*r4*r7*r6)/pow(rho,7.0/3.0));
      r14 = 1/pow(6,2.0/3.0);
      r15 = 1/pow(pi,4.0/3.0);
      r16 = 1/pow(rho,11.0/3.0);
      r17 = 1/pow(rho,8.0/3.0);
      r18 = -((1.0/3.0*r14*a*b*r15*r16*gamma)/sqrt((r14*pow(b,2.0)*r15*r17
                                                    *gamma)/4.0+1));
      r19 = 1/pow(rho,19.0/3.0);
      r20 = exp(-25.0*r14*r15*r17*gamma);
      r21 = -(d*r20)+c;
      r22 = (r14*r15*r17*gamma*r21)/4.0+r8+1;
      dpw_dr = -((r0*r1*a_x*r2*(r18-(1.0/4500.0*r3*r9*r19*r10)+r13)*r22)/
      pow(r11,2.0))+4.0/3.0*r0*r1*a_x*pow(rho,1.0/3.0)*r12*r22+r0*r1*a_x*
       r2*r12*(-(2.0/3.0*r14*r15*r16*gamma*r21)-(25.0/9.0*r3*d*r9*r19*r10*
       r20)+r18+r13);
      double g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18;
      g0 = 1/pow(3,1.0/3.0);
      g1 = pow(6,1.0/3.0);
      g2 = pow(rho,4.0/3.0);
      g3 = 1/g1;
      g4 = 1/pow(pi,2.0/3.0);
      g5 = 1/g2;
      g6 = sqrt(gamma);
      g7 = asinh((g3*b*g4*g5*g6)/2.0);
      g8 = (g3*a*g4*g5*g7*g6)/2.0;
      g9 = 1/pow(pi,8.0/3.0);
      g10 = 1/pow(rho,16.0/3.0);
      g11 = (g3*g9*g10*pow(gamma,2.0))/24000.0+g8+1;
      g12 = (g3*a*g4*g5*g7)/4.0/g6;
      g13 = 1/pow(6,2.0/3.0);
      g14 = 1/pow(pi,4.0/3.0);
      g15 = 1/pow(rho,8.0/3.0);
      g16 = (g13*a*b*g14*g15)/8.0/sqrt((g13*pow(b,2.0)*g14*g15*gamma)/4.0+
                                       1);
      g17 = exp(-25.0*g13*g14*g15*gamma);
      g18 = -(d*g17)+c;
      dpw_dg = -((g0*g1*a_x*g2*(g16+(g3*g9*g10*gamma)/12000.0+g12)*((g13*
       g14*g15*gamma*g18)/4.0+g8+1))/pow(g11,2.0))+(g0*g1*a_x*g2*((g13*g14
       *g15*g18)/4.0+25.0/24.0*g3*d*g9*g10*gamma*g17+g16+g12))/g11;
    }
}

void
PW91XFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double mpw, dmpw_dr, dmpw_dg;

  // alpha
  spin_contrib(id.a, mpw, dmpw_dr, dmpw_dg);
  od.energy = mpw;
  od.df_drho_a = dmpw_dr;
  od.df_dgamma_aa = dmpw_dg;

  // beta
  if (spin_polarized_) {
      spin_contrib(id.b, mpw, dmpw_dr, dmpw_dg);
      od.energy += mpw;
      od.df_drho_b = dmpw_dr;
      od.df_dgamma_bb = dmpw_dg;
    }
  else {
      od.energy += mpw;
      od.df_drho_b = od.df_drho_a;
      //Previously, we used just this:
      //od.df_dgamma_bb = od.df_dgamma_aa;
      //To be consistent with the other functionals, we use this:
      od.df_dgamma_aa *= 0.5;
      od.df_dgamma_bb = od.df_dgamma_aa;
      od.df_dgamma_ab = 2.*od.df_dgamma_aa;
    }
}

/////////////////////////////////////////////////////////////////////////////
// PBEXFunctional

static ClassDesc PBEXFunctional_cd(
  typeid(PBEXFunctional),"PBEXFunctional",1,"public DenFunctional",
  0, create<PBEXFunctional>, create<PBEXFunctional>);

PBEXFunctional::PBEXFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(mu);
  s.get(kappa);
}

PBEXFunctional::PBEXFunctional()
{
  init_constants();
}

PBEXFunctional::PBEXFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  mu = keyval->doublevalue("mu", KeyValValuedouble(mu));
  kappa = keyval->doublevalue("kappa", KeyValValuedouble(kappa));
  if (keyval->booleanvalue("revPBE")) {
      kappa = 1.245;
    }
}

PBEXFunctional::~PBEXFunctional()
{
}

void
PBEXFunctional::init_constants()
{
  // in paper:
  // mu = 0.21951;
  // in PBE.F
  mu = 0.2195149727645171;
  kappa = 0.804;
}

void
PBEXFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(mu);
  s.put(kappa);
}

int
PBEXFunctional::need_density_gradient()
{
  return 1;
}

void
PBEXFunctional::spin_contrib(const PointInputData::SpinData &i,
                             double &pbex,
                             double &dpbex_drhoa, double &dpbex_dgaa)
{
  double rhoa = i.rho;
  double rhoa_13 = i.rho_13;
  double rhoa_43 = rhoa*rhoa_13;
  double gaa = i.gamma;

  if (rhoa < MIN_DENSITY) {
      pbex = 0;
      if (compute_potential_) {
          dpbex_drhoa = 0.0;
          dpbex_dgaa = 0.0; // really -inf
        }
      return;
    }

  if (gaa < MIN_GAMMA) {
      pbex = - 0.93052573634910019 * rhoa_43;
      if (compute_potential_) {
          dpbex_drhoa = -(1.2407009817988002*rhoa_13);
          dpbex_dgaa = -((0.015312087450269402*mu)/rhoa_43);
        }
      return;
    }

  double rhoa_83 = rhoa_43*rhoa_43;

  pbex = -(0.93052573634910007*(-(kappa/((
      0.016455307846020562*gaa*mu)/kappa/rhoa_83+1.0))+kappa+1.0)*rhoa_43);

  if (compute_potential_) {
      double rhoa_73 = rhoa_43*rhoa;
      double r0;
      r0 = (0.016455307846020562*gaa*mu)/kappa/rhoa_83+1.0;
      dpbex_drhoa = -(1.2407009817988002*(-(kappa/r0)
       +kappa+1.0)*rhoa_13)+(0.040832233200718403*gaa*mu)/(r0*r0)/rhoa_73;
      dpbex_dgaa = -((0.015312087450269402*mu)/(r0*r0)/rhoa_43);
    }
}

void
PBEXFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double pbex, dpbex_dr, dpbex_dg;

  // alpha
  spin_contrib(id.a, pbex, dpbex_dr, dpbex_dg);
  od.energy = pbex;
  if (compute_potential_) {
      od.df_drho_a = dpbex_dr;
      od.df_dgamma_aa = dpbex_dg;
    }

  // beta
  if (spin_polarized_) {
      spin_contrib(id.b, pbex, dpbex_dr, dpbex_dg);
      od.energy += pbex;
      if (compute_potential_) {
          od.df_drho_b = dpbex_dr;
          od.df_dgamma_bb = dpbex_dg;
        }
    }
  else {
      od.energy += pbex;
      if (compute_potential_) {
          od.df_drho_b = od.df_drho_a;
          od.df_dgamma_bb = od.df_dgamma_aa;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// mPW91XFunctional

static ClassDesc mPW91XFunctional_cd(
  typeid(mPW91XFunctional),"mPW91XFunctional",1,"public DenFunctional",
  0, create<mPW91XFunctional>, create<mPW91XFunctional>);

mPW91XFunctional::mPW91XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  init_constants(mPW91);
  s.get(b);
  s.get(beta);
  s.get(c);
  s.get(d);
  s.get(x_d_coef);
}

mPW91XFunctional::mPW91XFunctional()
{
  init_constants(mPW91);
}

mPW91XFunctional::mPW91XFunctional(mPW91XFunctional::Func f)
{
  init_constants(f);
}

mPW91XFunctional::mPW91XFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  if (keyval->exists("constants")) {
      std::string t = keyval->stringvalue("constants");
      if (t == "B88") {
          init_constants(B88);
        }
      else if (t == "PW91") {
          init_constants(PW91);
        }
      else if (t == "mPW91") {
          init_constants(mPW91);
        }
      else {
          ExEnv::outn() << "mPW91XFunctional: bad \"constants\": " << t << std::endl;
          abort();
        }
    }
  else {
      init_constants(mPW91);
    }
  b = keyval->doublevalue("b", KeyValValuedouble(b));
  beta = keyval->doublevalue("beta", KeyValValuedouble(beta));
  c = keyval->doublevalue("c", KeyValValuedouble(c));
  d = keyval->doublevalue("d", KeyValValuedouble(d));
  x_d_coef = keyval->doublevalue("x_d_coef", KeyValValuedouble(x_d_coef));
}

mPW91XFunctional::~mPW91XFunctional()
{
}

void
mPW91XFunctional::init_constants(Func f)
{
  a_x = -1.5*pow(3./(4.*M_PI), 1./3.);
  if (f == B88) {
      b = 0.0042;
      beta = 0.0042;
      c = 0.;
      d = 1.;
      x_d_coef = 0.;
    }
  else if (f == PW91) {
      b = 0.0042;
      beta = 5.*pow(36.*M_PI,-5./3.);
      c = 1.6455;
      d = 4.;
      x_d_coef = 1.e-6;
    }
  else {
      b = 0.0046;
      beta = 5.*pow(36.*M_PI,-5./3.);
      c = 1.6455;
      d = 3.73;
      x_d_coef = 1.e-6;
    }
}

void
mPW91XFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(b);
  s.put(beta);
  s.put(c);
  s.put(d);
  s.put(x_d_coef);
}

int
mPW91XFunctional::need_density_gradient()
{
  return 1;
}

void
mPW91XFunctional::spin_contrib(const PointInputData::SpinData &i,
                               double &mpw, double &dmpw_dr, double &dmpw_dg)
{
  // this has been generated by macsyma
  double rho = i.rho;
  double gamma = i.gamma;

  if (rho < MIN_DENSITY) {
      mpw = 0.;
      dmpw_dr = 0.;
      dmpw_dg = 0.; // division by zero
      return;
    }
  if (gamma < MIN_GAMMA) {
      double rho_43 = rho * i.rho_13; // rho^(4/3)
      // 2^(1/3) * a_x * rho^(4/3)
      mpw = a_x * rho_43;
      // (4/3) a_x rho^(1/3)
      dmpw_dr = (4./3.) * a_x * i.rho_13;
      dmpw_dg = -beta/rho_43;
      return;
    }

  double t0,t1,t2,t3,t4,t5;
  t0 = pow(rho,4.0/3.0);
  t1 = 1/t0;
  t2 = sqrt(gamma);
  t3 = 1/pow(rho,4.0/3.0*d);
  t4 = pow(gamma,d/2.0);
  t5 = 1/pow(rho,8.0/3.0);
  mpw = t0*(-((-(((-beta+b)*t5*gamma)*exp(-c*t5*gamma))-(t3*x_d_coef*t4
     )+b*t5*gamma)/(-((t3*x_d_coef*t4)/a_x)+6*b*t1*asinh(t1*t2)*t2+1))+
     a_x);
  if (compute_potential_) {
      double r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14;
      r0 = pow(rho,4.0/3.0);
      r1 = 1/r0;
      r2 = sqrt(gamma);
      r3 = asinh(r1*r2);
      r4 = 1/a_x;
      r5 = 1/pow(rho,4.0/3.0*d);
      r6 = pow(gamma,d/2.0);
      r7 = -(r4*r5*x_d_coef*r6)+6*b*r1*r3*r2+1;
      r8 = 1/r7;
      r9 = 1/pow(rho,8.0/3.0);
      r10 = -beta+b;
      r11 = exp(-c*r9*gamma);
      r12 = -(r10*r9*gamma*r11)-(r5*x_d_coef*r6)+b*r9*gamma;
      r13 = 1/pow(rho,11.0/3.0);
      r14 = pow(rho,-(4.0/3.0*d)-1);
      dmpw_dr = r0*(-(r8*(-((8.0/3.0*r10*c*pow(gamma,2.0)*r11)/pow(rho,
       19.0/3.0))+8.0/3.0*r10*r13*gamma*r11+4.0/3.0*d*r14*x_d_coef*r6-(8.0
       /3.0*b*r13*gamma)))+((4.0/3.0*r4*d*r14*x_d_coef*r6-((8*b*r13*gamma)
       /sqrt(r9*gamma+1))-((8*b*r3*r2)/pow(rho,7.0/3.0)))*r12)/pow(r7,2.0)
       )+4.0/3.0*pow(rho,1.0/3.0)*(-(r8*r12)+a_x);
      double g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12;
      g0 = pow(rho,4.0/3.0);
      g1 = 1/g0;
      g2 = sqrt(gamma);
      g3 = asinh(g1*g2);
      g4 = 1/a_x;
      g5 = 1/pow(rho,4.0/3.0*d);
      g6 = d/2.0;
      g7 = pow(gamma,g6);
      g8 = -(g4*g5*x_d_coef*g7)+6*b*g1*g3*g2+1;
      g9 = 1/pow(rho,8.0/3.0);
      g10 = pow(gamma,g6-1);
      g11 = -beta+b;
      g12 = exp(-c*g9*gamma);
      dmpw_dg = g0*(((-(1.0/2.0*g4*d*g5*x_d_coef*g10)+(3*b*g9)/sqrt(g9*
       gamma+1)+(3*b*g1*g3)/g2)*(-(g11*g9*gamma*g12)-(g5*x_d_coef*g7)+b*g9
       *gamma))/pow(g8,2.0)-(((g11*c*gamma*g12)/pow(rho,16.0/3.0)-(g11*g9*
       g12)-(1.0/2.0*d*g5*x_d_coef*g10)+b*g9)/g8));
    }
}

void
mPW91XFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double mpw, dmpw_dr, dmpw_dg;

  // alpha
  spin_contrib(id.a, mpw, dmpw_dr, dmpw_dg);
  od.energy = mpw;
  od.df_drho_a = dmpw_dr;
  od.df_dgamma_aa = dmpw_dg;

  // beta
  if (spin_polarized_) {
      spin_contrib(id.b, mpw, dmpw_dr, dmpw_dg);
      od.energy += mpw;
      od.df_drho_b = dmpw_dr;
      od.df_dgamma_bb = dmpw_dg;
    }
  else {
      od.energy += mpw;
      od.df_drho_b = od.df_drho_a;
      od.df_dgamma_bb = od.df_dgamma_aa;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Wang (PW86) Exchange Functional
// J. P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800, 1986. 
// 
// Coded by Matt Leininger

static ClassDesc PW86XFunctional_cd(
  typeid(PW86XFunctional),"PW86XFunctional",1,"public DenFunctional",
  0, create<PW86XFunctional>, create<PW86XFunctional>);

PW86XFunctional::PW86XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(m_);
  s.get(a_);
  s.get(b_);
  s.get(c_);
}

PW86XFunctional::PW86XFunctional()
{
  init_constants();
}

PW86XFunctional::PW86XFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  m_ = keyval->doublevalue("m", KeyValValuedouble(m_));
  a_ = keyval->doublevalue("a", KeyValValuedouble(a_));
  b_ = keyval->doublevalue("b", KeyValValuedouble(b_));
  c_ = keyval->doublevalue("c", KeyValValuedouble(c_));
}

PW86XFunctional::~PW86XFunctional()
{
}

void
PW86XFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(m_);
  s.put(a_);
  s.put(b_);
  s.put(c_);
}

void
PW86XFunctional::init_constants()
{
  m_ = 1./15.;
  a_ = 0.0864/m_;
  b_ = 14.;
  c_ = 0.2;
}

int
PW86XFunctional::need_density_gradient()
{
  return 1;
}

void
 PW86XFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double rhoa = 2. * id.a.rho;
  double k_fa = pow( (3.*M_PI*M_PI*rhoa), (1./3.) );
  double rhoa43 = pow(rhoa, (4./3.));
  double rhoa13 = pow(rhoa, (1./3.));
  double Ax = -3./4.*pow( (3./M_PI), (1./3.) );
  double gamma_aa = 2. * sqrt(id.a.gamma);
  double sa;
  if (rhoa < MIN_DENSITY) sa = 0.;
  else sa = gamma_aa/(2. * k_fa * rhoa);
  double sa2 = sa*sa;
  double sa3 = sa2*sa;
  double sa4 = sa2*sa2;
  double sa5 = sa4*sa;
  double sa6 = sa5*sa;
  double F1a = 1. + a_*sa2 + b_*sa4 + c_*sa6;
  double Fxa = pow(F1a, m_);
  double fpw86xa = Ax * rhoa43 * Fxa;
  double ex = 0.5 * fpw86xa;

  if (compute_potential_) {
    double dsa_drhoa;
    if (rhoa < MIN_DENSITY) dsa_drhoa = 0.;
    else dsa_drhoa = -4.*2.*sa/(3.*rhoa);
    double dFxa_drhoa = m_ * pow(F1a, (m_-1.)) * ( 2.*a_*sa + 4.*b_*sa3 + 6.*c_*sa5) * dsa_drhoa;
    double dfpw86xa_drhoa = Ax * ( 2.*4./3.*rhoa13*Fxa + rhoa43 * dFxa_drhoa );
    od.df_drho_a = 0.5 * dfpw86xa_drhoa; 

    double sa_dsa_dgamma_aa;
    if (rhoa < MIN_DENSITY) sa_dsa_dgamma_aa = 0.;
    else {
        sa_dsa_dgamma_aa = 2./pow( (2.*k_fa*rhoa), 2.);
      }
    double dFxa_dgamma_aa = m_ * pow(F1a, (m_-1.)) *
                    ( 2.*a_ + 4.*b_*sa2 + 6.*c_*sa4 ) * sa_dsa_dgamma_aa ;
    double dfpw86xa_dgamma_aa = Ax * rhoa43 * dFxa_dgamma_aa;
    od.df_dgamma_aa = 0.5 * dfpw86xa_dgamma_aa;

      od.df_drho_b = od.df_drho_a;
      od.df_dgamma_bb = od.df_dgamma_aa;
      od.df_dgamma_ab = 0.;
    }
  
    if (spin_polarized_) {
      double rhob = 2. * id.b.rho;
      double k_fb = pow( (3.*M_PI*M_PI*rhob), (1./3.) );
      double rhob43 = pow(rhob, (4./3.));
      double rhob13 = pow(rhob, (1./3.));
      double gamma_bb = 2.*sqrt(id.b.gamma);
      double sb;
      if (rhob < MIN_DENSITY) sb = 0.;
      else sb = gamma_bb/(2.*k_fb*rhob);
      double sb2 = sb*sb;
      double sb3 = sb2*sb;
      double sb4 = sb3*sb;
      double sb5 = sb4*sb;
      double sb6 = sb5*sb;
      double F1b = 1. + a_*sb2 + b_*sb4 + c_*sb6;
      double Fxb = pow(F1b, m_);
      double fpw86xb = Ax * rhob43 * Fxb;
      ex += 0.5 * fpw86xb;
      
      if (compute_potential_) {
          double dsb_drhob;
          if (rhob < MIN_DENSITY) dsb_drhob = 0.;
          else dsb_drhob = -4.*2.*sb/(3.*rhob);
          double dFxb_drhob = m_ * pow(F1b, (m_-1.))
                             * (2.*a_*sb + 4.*b_*sb3 + 6.*c_*sb5) * dsb_drhob;
          double dfpw86xb = Ax * ( 2.*4./3.*rhob13*Fxb + rhob43 * dFxb_drhob );
          od.df_drho_b = 0.5 * dfpw86xb;

          double sb_dsb_dgamma_bb;
          if (rhob < MIN_DENSITY) sb_dsb_dgamma_bb = 0.;
          else {
              sb_dsb_dgamma_bb = 2./pow( (2.*k_fb*rhob), 2.);
            }
          double dFxb_dgamma_bb = m_ * pow(F1b, (m_-1.)) *
                    ( 2.*a_ + 4.*b_*sb2 + 6.*c_*sb4) * sb_dsb_dgamma_bb;
          double dfpw86xb_dgamma_bb = Ax * rhob43 * dFxb_dgamma_bb;
          od.df_dgamma_bb = 0.5 * dfpw86xb_dgamma_bb;
       }
    }
  else ex += ex;

  od.energy = ex;
}

/////////////////////////////////////////////////////////////////////////////
// Gill 1996 (G96) Exchange Functional
// P. M. W. Gill, Mol. Phys. 89, 433, 1996.
// 
// Coded by Matt Leininger

static ClassDesc G96XFunctional_cd(
  typeid(G96XFunctional),"G96XFunctional",1,"public DenFunctional",
  0, create<G96XFunctional>, create<G96XFunctional>);

G96XFunctional::G96XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
  s.get(b_);
}

G96XFunctional::G96XFunctional()
{
  init_constants();
}

G96XFunctional::G96XFunctional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
  init_constants();
  b_ = keyval->doublevalue("b", KeyValValuedouble(b_));
}

G96XFunctional::~G96XFunctional()
{
}

void
G96XFunctional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
  s.put(b_);
}

void
G96XFunctional::init_constants()
{
  b_ = 1./137.;
}

int
G96XFunctional::need_density_gradient()
{
  return 1;
}

void
G96XFunctional::point(const PointInputData &id,
                      PointOutputData &od)
{
  od.zero();

  double rhoa = id.a.rho;
  double rhoa43 = pow(rhoa, (4./3.));
  double rhoa13 = pow(rhoa, (1./3.));
  double gamma_aa = sqrt(id.a.gamma);
  double gamma_aa32 = pow(gamma_aa, 1.5);
  double alpha = -1.5 * pow( (3./(4.*M_PI)), (1./3.) );
  double gg96a;
  double fxg96a;
  if (rhoa < MIN_DENSITY) gg96a = fxg96a = 0.;
  else {
      gg96a = alpha - b_*gamma_aa32/(rhoa*rhoa);
      fxg96a = rhoa43 * gg96a;
    }
  double ex = fxg96a;

  if (compute_potential_) {

      double dfxg96a_drhoa;
      if (rhoa < MIN_DENSITY) dfxg96a_drhoa = 0.;
      else {
          double dgg96a_drhoa = 2.*b_*gamma_aa32/(rhoa*rhoa*rhoa);
          dfxg96a_drhoa = 4./3.*rhoa13*gg96a + rhoa43*dgg96a_drhoa;
        }
      od.df_drho_a = dfxg96a_drhoa;
      od.df_drho_b = od.df_drho_a;
      
      double dfxg96a_dgamma_aa;
      // The derivative of the G96X functional with respect to gamma_aa or bb
      // as implemented should go to infinity as gamma goes to zero.
      // However, the derivative gamma terms are eventually contracted with quantities
      // that have a sqrt(gamma) in the numerator and therefore the overall limit
      // is zero.
      if (gamma_aa < MIN_GAMMA || rhoa < MIN_DENSITY ) dfxg96a_dgamma_aa = 0.;
      else {
          double dgg96a_dgamma_aa = -3.*b_ / ( 4.*rhoa*rhoa*sqrt(gamma_aa) );
          dfxg96a_dgamma_aa = rhoa43 * dgg96a_dgamma_aa;
        }
      od.df_dgamma_aa = dfxg96a_dgamma_aa;
      od.df_dgamma_bb = od.df_dgamma_aa;
      od.df_dgamma_ab = 0.;
    }
  
  if (spin_polarized_) {
      double rhob = id.b.rho;
      double rhob43 = pow(rhob, (4./3.));
      double rhob13 = pow(rhob, (1./3.));
      double gamma_bb = sqrt(id.b.gamma);
      double gamma_bb32 = pow(gamma_bb, 1.5);
      double gg96b;
      double fxg96b;
      if (rhob < MIN_DENSITY) gg96b = fxg96b = 0.;
      else {
          gg96b = alpha - b_*gamma_bb32/(rhob*rhob);
          fxg96b = rhob43 * gg96b;
        }
      ex += fxg96b;
    
      if (compute_potential_) {
          double dfxg96b_drhob;
          if (rhob < MIN_DENSITY) dfxg96b_drhob = 0.;
          else {
              double dgg96b_drhob = 2.*b_*gamma_bb32/(rhob*rhob*rhob);
              dfxg96b_drhob = 4./3.*rhob13*gg96b + rhob43*dgg96b_drhob;
            }
          od.df_drho_b = dfxg96b_drhob;

          double dfxg96b_dgamma_bb;
          // See comment above with regard to correct limits.
          if (gamma_bb < MIN_GAMMA || rhob < MIN_DENSITY) dfxg96b_dgamma_bb=0.;
          else {
              double dgg96b_dgamma_bb = -3.*b_ / (4.*rhob*rhob*sqrt(gamma_bb));
              dfxg96b_dgamma_bb = rhob43 * dgg96b_dgamma_bb;
            }
          od.df_dgamma_bb = dfxg96b_dgamma_bb;
        }
      }
  else ex += ex;
  
  
  od.energy = ex;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
