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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/dft/functional.h>

///////////////////////////////////////////////////////////////////////////
// PointInputData

void
PointInputData::compute_derived(int spin_polarized)
{
  a.rho_13 = pow(a.rho, 1.0/3.0);
  if (spin_polarized) {
      b.rho_13 = pow(b.rho, 1.0/3.0);
    }
  else {
      b = a;
      gamma_ab = a.gamma;
    }
}

///////////////////////////////////////////////////////////////////////////
// DenFunctional

#define CLASSNAME DenFunctional
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
DenFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

DenFunctional::DenFunctional(StateIn& s):
  SavableState(s)
{
  s.get(a0_);
}

DenFunctional::DenFunctional()
{
  a0_ = 0;
  spin_polarized_ = 0;
  compute_potential_ = 0;
}

DenFunctional::DenFunctional(const RefKeyVal& keyval)
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
}

int
DenFunctional::need_density_gradient()
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
                        const RefGaussianBasisSet &basis,
                        const double *dmat_a, const double *dmat_b,
                        const double *bs, const double *bsg,
                        const double *bsh)
{
  int need_gamma_terms = need_density_gradient();
  point(id, od);
  memset(grad_f, 0, sizeof(double)*basis->molecule()->natom()*3);
#if 0
  cout << scprintf("gradient: rho_a= %12.8f rho_b= %12.8f need_gamma = %d",
                   id.a.rho, id.b.rho, need_gamma_terms) << endl;
  cout << scprintf("  gamma_aa= %12.8f gamma_bb= %12.8f gamma_ab= % 12.8f",
                   id.a.gamma, id.b.gamma, id.gamma_ab) << endl;
  cout << scprintf("  df_drho_a= % 12.8f df_drho_b= % 12.8f",
                   od.df_drho_a, od.df_drho_b) << endl;
  cout << scprintf("  df_dg_aa= % 12.8f df_dg_bb= % 12.8f df_dg_ab= % 12.8f",
                   od.df_dgamma_aa,od.df_dgamma_bb,od.df_dgamma_ab) << endl;
#endif
  int nbasis = basis->nbasis();
  int jk=0;
  double drhoa = od.df_drho_a;
  double drhob = od.df_drho_b;
  for (int nu=0; nu < nbasis; nu++) {
      int nuatom
          = basis->shell_to_center(basis->function_to_shell(nu));
      double dfa_phi_nu = drhoa * bs[nu];
      double dfb_phi_nu = drhob * bs[nu];
      int nush = basis->function_to_shell(nu);
      for (int mu=0; mu<nbasis; mu++) {
          int muatom
              = basis->shell_to_center(basis->function_to_shell(mu));
          if (muatom!=acenter) {
              int numu = (nu>mu?((nu*(nu+1))/2+mu):((mu*(mu+1))/2+nu));
              double rho_a = dmat_a[numu];
              double rho_b = dmat_b[numu];
              int ixyz;
              for (ixyz=0; ixyz<3; ixyz++) {
                  double contrib = -2.0*bsg[mu*3+ixyz]
                                 * (rho_a*dfa_phi_nu + rho_b*dfb_phi_nu);
#define hoff(i,j) ((j)<(i)?((i)*((i)+1))/2+(j):((j)*((j)+1))/2+(i))
                  // gamma_aa contrib
                  if (need_gamma_terms) {
                      contrib += 4.0 * od.df_dgamma_aa * rho_a
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
                      contrib += 2.0 * od.df_dgamma_ab * rho_a
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
                      contrib += 2.0 * od.df_dgamma_ab * rho_b
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
                      contrib += 4.0 * od.df_dgamma_bb * rho_b
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

void
DenFunctional::do_fd_point(PointInputData&id,
                           double&in,double&out,
                           double lower_bound, double upper_bound)
{
  double delta = 0.0000001;
  PointOutputData tod;
  double insave = in;

  point(id,tod);
  double outsave = tod.energy;

  int spin_polarized_save = spin_polarized_;
  set_spin_polarized(1);

  if (insave-delta>=lower_bound && insave+delta<=upper_bound) {
      in = insave+delta;
      id.compute_derived(1);
      point(id,tod);
      double plus = tod.energy;

      in = insave-delta;
      id.compute_derived(1);
      point(id,tod);
      double minu = tod.energy;
      out = 0.5*(plus-minu)/delta;
    }
  else if (insave+2*delta<=upper_bound) {
      in = insave+delta;
      id.compute_derived(1);
      point(id,tod);
      double plus = tod.energy;

      in = insave+2*delta;
      id.compute_derived(1);
      point(id,tod);
      double plus2 = tod.energy;
      out = 0.5*(4.0*plus-plus2-3.0*outsave)/delta;
    }
  else if (insave-2*delta>=lower_bound) {
      in = insave-delta;
      id.compute_derived(1);
      point(id,tod);
      double minu = tod.energy;

      in = insave-2*delta;
      id.compute_derived(1);
      point(id,tod);
      double minu2 = tod.energy;
      out = -0.5*(4.0*minu-minu2-3.0*outsave)/delta;
    }
  else {
      // the derivative is not well defined for this case
      out = 0.0;
    }
  in = insave;
  id.compute_derived(1);

  set_spin_polarized(spin_polarized_save);
}

void
DenFunctional::fd_point(const PointInputData&id, PointOutputData&od)
{
  PointInputData tid(id);
  double plus, minu;

  // fill in the energy at the initial density values
  point(id,od);

  cout << scprintf("ra= %6.4f rb= %6.4f gaa= %6.4f gbb= %6.4f gab= % 6.4f",
                   id.a.rho, id.b.rho, id.a.gamma, id.b.gamma, id.gamma_ab)
       << endl;

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

static void
check(const char *name, double fd, double an)
{
  double err = fabs(fd - an);
  cout << scprintf("%20s: fd = % 12.8f an = % 12.8f", name, fd, an)
       << endl;
  if ((fabs(an) > 0.03 && err/fabs(an) > 0.03)
      || ((fabs(an) <= 0.03) && err > 0.03)
      || isnan(an)) {
      cout << scprintf("Error: %20s: fd = % 12.8f an = % 12.8f", name, fd, an)
           << endl;
    }
}

void
DenFunctional::test(const PointInputData &id)
{
  PointOutputData fd_od;
  fd_point(id,fd_od);
  PointOutputData an_od;
  point(id,an_od);
  check("df_drho_a", fd_od.df_drho_a, an_od.df_drho_a);
  check("df_drho_b", fd_od.df_drho_b, an_od.df_drho_b);
  check("df_dgamma_aa", fd_od.df_dgamma_aa, an_od.df_dgamma_aa);
  check("df_dgamma_ab", fd_od.df_dgamma_ab, an_od.df_dgamma_ab);
  check("df_dgamma_bb", fd_od.df_dgamma_bb, an_od.df_dgamma_bb);
}

void
DenFunctional::test()
{
  int i, j, k, l, m;
  set_compute_potential(1);
  set_spin_polarized(0);
  SCVector3 r = 0.0;
  PointInputData id(r);

  for (i=0; i<6; i++) id.a.hes_rho[i] = 0.0;
  id.a.lap_rho = 0.0;

  // del rho should not be used by any of the functionals
  for (i=0; i<3; i++) id.a.del_rho[i] = id.b.del_rho[i] = 0.0;

  double testrho[] = { 0.001, 0.5, -1 };
  double testgamma[] = { 0.0001, 0.001, 0.5, -1 };
  double testgammaab[] = { -0.5, 0.0, 0.5, -1 };

  cout << "Testing with rho_a == rho_b" << endl;
  for (i=0; testrho[i] != -1.0; i++) {
      id.a.rho=testrho[i];
      for (j=0; testgamma[j] != -1.0; j++) {
          id.a.gamma = testgamma[j];
          id.compute_derived(0);
          test(id);
        }
    }

  set_spin_polarized(1);
  cout << "Testing with rho_a != rho_b" << endl;
  for (i=0; testrho[i] != -1.0; i++) {
      id.a.rho=testrho[i];
      for (j=0; testrho[j] != -1.0; j++) {
          id.b.rho=testrho[j];
          for (k=0; testgamma[k] != -1.0; k++) {
              id.a.gamma = testgamma[k];
              double sqrt_gamma_a = sqrt(id.a.gamma);
              for (l=0; testgamma[l] != -1.0; l++) {
                  id.b.gamma = testgamma[l];
                  double sqrt_gamma_b = sqrt(id.b.gamma);
                  for (m=0; testgammaab[m] != -1.0; m++) {
                      // constrain gamma_ab to values allowed by the
                      // current gamma_a and gamma_b
                      id.gamma_ab = testgammaab[m];
                      if (id.gamma_ab > sqrt_gamma_a*sqrt_gamma_b) {
                          id.gamma_ab = sqrt_gamma_a*sqrt_gamma_b;
                        }
                      if (id.gamma_ab < -0.5*(id.a.gamma+id.b.gamma)) {
                          id.gamma_ab = -0.5*(id.a.gamma+id.b.gamma);
                        }
                      if (id.gamma_ab < -sqrt_gamma_a*sqrt_gamma_b) {
                          id.gamma_ab = -sqrt_gamma_a*sqrt_gamma_b;
                        }
                      id.compute_derived(1);
                      test(id);
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// NElFunctional

#define CLASSNAME NElFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NElFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

NElFunctional::NElFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

NElFunctional::NElFunctional(const RefKeyVal& keyval):
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

#define CLASSNAME SumDenFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SumDenFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

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
      funcs_ = new RefDenFunctional[n_];
      for (int i=0; i < n_; i++)
          funcs_[i].restore_state(s);
    }
}

SumDenFunctional::SumDenFunctional() :
  n_(0),
  funcs_(0),
  coefs_(0)
{
}

SumDenFunctional::SumDenFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval),
  n_(0),
  funcs_(0),
  coefs_(0)
{
  int ncoef = keyval->count("coefs");
  int nfunc = keyval->count("funcs");
  if (ncoef != nfunc && ncoef != 0) {
      cout << "SumDenFunctional: number of coefs and funcs differ" << endl;
      abort();
    }
  
  n_ = nfunc;
  coefs_ = new double[n_];
  funcs_ = new RefDenFunctional[n_];
  for (int i=0; i < n_; i++) {
      if (ncoef)
          coefs_[i] = keyval->doublevalue("coefs", i);
      else
          coefs_[i] = 1.0;
      funcs_[i] = keyval->describedclassvalue("funcs", i);
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
  s.put(n_);
  if (n_) {
      s.put(coefs_, n_);
      for (int i=0; i < n_; i++) 
          funcs_[i].save_state(s);
    }
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
SumDenFunctional::print(ostream& o) const
{
  o << node0
    << indent << "Sum of Functionals:" << endl;
  o << incindent;
  for (int i=0; i<n_; i++) {
      o << node0 << indent << scprintf("%+18.16f",coefs_[i]) << endl;
      o << incindent;
      funcs_[i]->print(o);
      o << decindent;
    }
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////
// LSDACFunctional: All local correlation functional inherit from this class.
// Coded by Matt Leininger
#define CLASSNAME LSDACFunctional
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
LSDACFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

LSDACFunctional::LSDACFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

LSDACFunctional::LSDACFunctional()
{
}

LSDACFunctional::LSDACFunctional(const RefKeyVal& keyval):
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

#define CLASSNAME SlaterXFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SlaterXFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

SlaterXFunctional::SlaterXFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

SlaterXFunctional::SlaterXFunctional()
{
}

SlaterXFunctional::SlaterXFunctional(const RefKeyVal& keyval):
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
#define CLASSNAME PW92LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW92LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW92LCFunctional::PW92LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

PW92LCFunctional::PW92LCFunctional()
{
}

PW92LCFunctional::PW92LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
}

PW92LCFunctional::~PW92LCFunctional()
{
}

void
PW92LCFunctional::save_data_state(StateOut& s)
{
  cout << "PW92LCFunctional: cannot save state" << endl;
  abort();
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
#define CLASSNAME PZ81LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PZ81LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PZ81LCFunctional::PZ81LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

PZ81LCFunctional::PZ81LCFunctional()
{
}

PZ81LCFunctional::PZ81LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
}

PZ81LCFunctional::~PZ81LCFunctional()
{
}

void
PZ81LCFunctional::save_data_state(StateOut& s)
{
  cout << "PZ81LCFunctional: cannot save state" << endl;
  abort();
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
// VWN1LCFunctional
// Coded by Matt Leininger

#define CLASSNAME VWN1LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN1LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN1LCFunctional::VWN1LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

VWN1LCFunctional::VWN1LCFunctional()
{
  Ap_ = 0.0310907;
  x0p_ = -0.10498;
  bp_ = 3.72744;
  cp_ = 12.9352;
  Af_ = 0.01554535;
  x0f_ = -0.32500;
  bf_ = 7.06042;
  cf_ = 18.0578;
}

VWN1LCFunctional::VWN1LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
    Ap_ = keyval->doublevalue("Ap", KeyValValuedouble(0.0310907));
    x0p_ = keyval->doublevalue("x0p", KeyValValuedouble(-0.10498));
    bp_ = keyval->doublevalue("bp", KeyValValuedouble(3.72744));
    cp_ = keyval->doublevalue("cp", KeyValValuedouble(12.9352));
    Af_ = keyval->doublevalue("Af", KeyValValuedouble(0.01554535));
    x0f_ = keyval->doublevalue("x0f", KeyValValuedouble(-0.32500));
    bf_ = keyval->doublevalue("bf", KeyValValuedouble(7.06042));
    cf_ = keyval->doublevalue("cf", KeyValValuedouble(18.0578));

    int vwn1rpa = keyval->intvalue("vwn1rpa", KeyValValueint(0));
    if (vwn1rpa) {
        x0p_ = -0.409286;
        bp_ = 13.0720;
        cp_ = 42.7198;
        x0f_ = -0.743294;
        bf_ = 20.1231;
        cf_ = 101.578;
      }
}

VWN1LCFunctional::~VWN1LCFunctional()
{
}

void
VWN1LCFunctional::save_data_state(StateOut& s)
{
  cout << "VWN1LCFunctional: cannot save state" << endl;
  abort();
}

double
VWN1LCFunctional::F(double x, double A, double x0, double b, double c)
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
VWN1LCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
      double zeta2 = zeta*zeta;
      double zeta3 = zeta2*zeta;
      double zeta4 = zeta2*zeta2;
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
#define CLASSNAME VWN2LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN2LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN2LCFunctional::VWN2LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

VWN2LCFunctional::VWN2LCFunctional()
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

  x0p_rpa_ = -0.409286;
  bp_rpa_  = 13.0720;
  cp_rpa_  = 42.7198;
  x0f_rpa_ = -0.743294;
  bf_rpa_  = 20.1231;
  cf_rpa_  = 101.578;

  x0_alpha_mc_ = -0.00475840;
  b_alpha_mc_  = 1.13107;
  c_alpha_mc_  = 13.0045;

  x0_alpha_rpa_ = -0.228344;
  b_alpha_rpa_  = 1.06835;
  c_alpha_rpa_  = 11.4813;
  
}

VWN2LCFunctional::VWN2LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
    Ap_  = keyval->doublevalue("Ap", KeyValValuedouble(0.0310907));
    Af_  = keyval->doublevalue("Af", KeyValValuedouble(0.01554535));
    A_alpha_  = keyval->doublevalue("A_alpha", KeyValValuedouble(-1./(6.*M_PI*M_PI)));
 
    x0p_mc_  = keyval->doublevalue("x0p_mc", KeyValValuedouble(-0.10498));
    bp_mc_   = keyval->doublevalue("bp_mc", KeyValValuedouble(3.72744));
    cp_mc_   = keyval->doublevalue("cp_mc", KeyValValuedouble(12.9352));
    x0f_mc_  = keyval->doublevalue("x0f_mc", KeyValValuedouble(-0.32500));
    bf_mc_   = keyval->doublevalue("bf_mc", KeyValValuedouble(7.06042));
    cf_mc_   = keyval->doublevalue("cf_mc", KeyValValuedouble(18.0578));

    x0p_rpa_ = keyval->doublevalue("x0p_rpa", KeyValValuedouble(-0.409286));
    bp_rpa_  = keyval->doublevalue("bp_rpa", KeyValValuedouble(13.0720));
    cp_rpa_  = keyval->doublevalue("cp_rpa", KeyValValuedouble(42.7198));
    x0f_rpa_ = keyval->doublevalue("x0f_rpa", KeyValValuedouble(-0.743294));
    bf_rpa_  = keyval->doublevalue("bf_rpa", KeyValValuedouble(20.1231));
    cf_rpa_  = keyval->doublevalue("cf_rpa", KeyValValuedouble(101.578));

    x0_alpha_mc_ = keyval->doublevalue("x0_alpha_mc", KeyValValuedouble(-0.00475840));
    b_alpha_mc_  = keyval->doublevalue("b_alpha_mc", KeyValValuedouble(1.13107));
    c_alpha_mc_  = keyval->doublevalue("c_alpha_mc", KeyValValuedouble(13.0045));

    x0_alpha_rpa_ = keyval->doublevalue("x0_alpha_rpa", KeyValValuedouble(-0.228344));
    b_alpha_rpa_  = keyval->doublevalue("b_alpha_rpa", KeyValValuedouble(1.06835));
    c_alpha_rpa_  = keyval->doublevalue("c_alpha_rpa", KeyValValuedouble(11.4813));
    
}

VWN2LCFunctional::~VWN2LCFunctional()
{
}

void
VWN2LCFunctional::save_data_state(StateOut& s)
{
  cout << "VWN2LCFunctional: cannot save state" << endl;
  abort();
}

double
VWN2LCFunctional::F(double x, double A, double x0, double b, double c)
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
VWN2LCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
  double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
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
  double ec = epc_rpa + delta_ec;

  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      double zeta3 = zeta2*zeta;
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
      
      double dec_dzeta = ddeltae_rpa_dzeta + fp * (ddelta_e_mc - ddelta_e_rpa);
                       
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN3LCFunctional
// Coded by Matt Leininger
#define CLASSNAME VWN3LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN3LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN3LCFunctional::VWN3LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

VWN3LCFunctional::VWN3LCFunctional()
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

  x0p_rpa_ = -0.409286;
  bp_rpa_  = 13.0720;
  cp_rpa_  = 42.7198;
  x0f_rpa_ = -0.743294;
  bf_rpa_  = 20.1231;
  cf_rpa_  = 101.578;

  x0_alpha_mc_ = -0.00475840;
  b_alpha_mc_  = 1.13107;
  c_alpha_mc_  = 13.0045;

  x0_alpha_rpa_ = -0.228344;
  b_alpha_rpa_  = 1.06835;
  c_alpha_rpa_  = 11.4813;

}

VWN3LCFunctional::VWN3LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
    Ap_  = keyval->doublevalue("Ap", KeyValValuedouble(0.0310907));
    Af_  = keyval->doublevalue("Af", KeyValValuedouble(0.01554535));
    A_alpha_  = keyval->doublevalue("A_alpha", KeyValValuedouble(-1./(6.*M_PI*M_PI)));

    x0p_mc_  = keyval->doublevalue("x0p_mc", KeyValValuedouble(-0.10498));
    bp_mc_   = keyval->doublevalue("bp_mc", KeyValValuedouble(3.72744));
    cp_mc_   = keyval->doublevalue("cp_mc", KeyValValuedouble(12.9352));
    x0f_mc_  = keyval->doublevalue("x0f_mc", KeyValValuedouble(-0.32500));
    bf_mc_   = keyval->doublevalue("bf_mc", KeyValValuedouble(7.06042));
    cf_mc_   = keyval->doublevalue("cf_mc", KeyValValuedouble(18.0578));

    x0p_rpa_ = keyval->doublevalue("x0p_rpa", KeyValValuedouble(-0.409286));
    bp_rpa_  = keyval->doublevalue("bp_rpa", KeyValValuedouble(13.0720));
    cp_rpa_  = keyval->doublevalue("cp_rpa", KeyValValuedouble(42.7198));
    x0f_rpa_ = keyval->doublevalue("x0f_rpa", KeyValValuedouble(-0.743294));
    bf_rpa_  = keyval->doublevalue("bf_rpa", KeyValValuedouble(20.1231));
    cf_rpa_  = keyval->doublevalue("cf_rpa", KeyValValuedouble(101.578));

    x0_alpha_mc_ = keyval->doublevalue("x0_alpha_mc", KeyValValuedouble(-0.00475840));
    b_alpha_mc_  = keyval->doublevalue("b_alpha_mc", KeyValValuedouble(1.13107));
    c_alpha_mc_  = keyval->doublevalue("c_alpha_mc", KeyValValuedouble(13.0045));

    x0_alpha_rpa_ = keyval->doublevalue("x0_alpha_rpa", KeyValValuedouble(-0.228344));
    b_alpha_rpa_  = keyval->doublevalue("b_alpha_rpa", KeyValValuedouble(1.06835));
    c_alpha_rpa_  = keyval->doublevalue("c_alpha_rpa", KeyValValuedouble(11.4813));

}

VWN3LCFunctional::~VWN3LCFunctional()
{
}

void
VWN3LCFunctional::save_data_state(StateOut& s)
{
  cout << "VWN3LCFunctional: cannot save state" << endl;
  abort();
}

double
VWN3LCFunctional::F(double x, double A, double x0, double b, double c)
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
VWN3LCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
  double epc_mc    = F(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
  double efc_mc    = F(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
  double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
  // RPA fitting parameters
  double epc_rpa    = F(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
  double efc_rpa    = F(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
  double alphac_rpa = F(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);

  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc_rpa - epc_rpa;
  double delta_e_mc  = efc_mc - epc_mc;
  double beta = fpp0 * delta_e_rpa / alphac_rpa - 1.;
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. + beta * zeta4);
  double delta_ec = delta_e_mc/delta_e_rpa * delta_erpa_rszeta;
  double ec = epc_rpa + delta_ec;
  // double ec = epc_mc + delta_ec;
  
  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      double zeta3 = zeta2*zeta;
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
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
      double ddelta_e_mc = defc_dr_s1_mc - depc_dr_s0_mc;
      // double dec_dr_s_mc = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
      //                   + dalphac_dr_s * f / fpp0 * (1 - zeta4);
      // double dec_dzeta_mc = 4.* zeta3 * f * (efc_mc - epc_mc - (alphac_mc/fpp0))
      //        + fp * (zeta4 * (efc_mc - epc_mc) + (1-zeta4)*(alphac_mc/fpp0));
      
      double dec_dr_s = depc_dr_s0_rpa + delta_erpa_rszeta/delta_e_rpa * ddelta_e_mc
        + delta_e_mc/delta_e_rpa * ddeltae_rpa_dr_s 
        - delta_erpa_rszeta*delta_e_mc/(delta_e_rpa*delta_e_rpa)* ddelta_e_rpa;
      double dec_dzeta = delta_e_mc / delta_e_rpa * ddeltae_rpa_dzeta;      
 
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN4LCFunctional
// Coded by Matt Leininger
#define CLASSNAME VWN4LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN4LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN4LCFunctional::VWN4LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

VWN4LCFunctional::VWN4LCFunctional()
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

  x0p_rpa_ = -0.409286;
  bp_rpa_  = 13.0720;
  cp_rpa_  = 42.7198;
  x0f_rpa_ = -0.743294;
  bf_rpa_  = 20.1231;
  cf_rpa_  = 101.578;

  x0_alpha_mc_ = -0.00475840;
  b_alpha_mc_  = 1.13107;
  c_alpha_mc_  = 13.0045;

  x0_alpha_rpa_ = -0.228344;
  b_alpha_rpa_  = 1.06835;
  c_alpha_rpa_  = 11.4813;

}

VWN4LCFunctional::VWN4LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
    Ap_  = keyval->doublevalue("Ap", KeyValValuedouble(0.0310907));
    Af_  = keyval->doublevalue("Af", KeyValValuedouble(0.01554535));
    A_alpha_  = keyval->doublevalue("A_alpha", KeyValValuedouble(-1./(6.*M_PI*M_PI)));

    x0p_mc_  = keyval->doublevalue("x0p_mc", KeyValValuedouble(-0.10498));
    bp_mc_   = keyval->doublevalue("bp_mc", KeyValValuedouble(3.72744));
    cp_mc_   = keyval->doublevalue("cp_mc", KeyValValuedouble(12.9352));
    x0f_mc_  = keyval->doublevalue("x0f_mc", KeyValValuedouble(-0.32500));
    bf_mc_   = keyval->doublevalue("bf_mc", KeyValValuedouble(7.06042));
    cf_mc_   = keyval->doublevalue("cf_mc", KeyValValuedouble(18.0578));

    x0p_rpa_ = keyval->doublevalue("x0p_rpa", KeyValValuedouble(-0.409286));
    bp_rpa_  = keyval->doublevalue("bp_rpa", KeyValValuedouble(13.0720));
    cp_rpa_  = keyval->doublevalue("cp_rpa", KeyValValuedouble(42.7198));
    x0f_rpa_ = keyval->doublevalue("x0f_rpa", KeyValValuedouble(-0.743294));
    bf_rpa_  = keyval->doublevalue("bf_rpa", KeyValValuedouble(20.1231));
    cf_rpa_  = keyval->doublevalue("cf_rpa", KeyValValuedouble(101.578));

    x0_alpha_mc_ = keyval->doublevalue("x0_alpha_mc", KeyValValuedouble(-0.00475840));
    b_alpha_mc_  = keyval->doublevalue("b_alpha_mc", KeyValValuedouble(1.13107));
    c_alpha_mc_  = keyval->doublevalue("c_alpha_mc", KeyValValuedouble(13.0045));

    x0_alpha_rpa_ = keyval->doublevalue("x0_alpha_rpa", KeyValValuedouble(-0.228344));
    b_alpha_rpa_  = keyval->doublevalue("b_alpha_rpa", KeyValValuedouble(1.06835));
    c_alpha_rpa_  = keyval->doublevalue("c_alpha_rpa", KeyValValuedouble(11.4813));

}

VWN4LCFunctional::~VWN4LCFunctional()
{
}

void
VWN4LCFunctional::save_data_state(StateOut& s)
{
  cout << "VWN4LCFunctional: cannot save state" << endl;
  abort();
}

double
VWN4LCFunctional::F(double x, double A, double x0, double b, double c)
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
VWN4LCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
  double alphac_mc = F(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
  // RPA fitting parameters
  double epc_rpa    = F(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
  double efc_rpa    = F(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
  double alphac_rpa = F(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_);

  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc_rpa - epc_rpa;
  double delta_e_mc  = efc_mc - epc_mc;
  double beta = fpp0 * delta_e_mc / alphac_rpa - 1.;
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. + beta * zeta4);
  double delta_ec = delta_e_mc/delta_e_rpa * delta_erpa_rszeta;
  // Is the epc used below suppose to be from the mc or rpa data?
  // We will use the mc data for now.
  double ec = epc_mc + delta_ec;

  od.energy = ec * rho;
  ec_local = ec;

  if (compute_potential_) {
      double zeta3 = zeta2*zeta;
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, Ap_, x0p_mc_, bp_mc_, cp_mc_);
      double defc_dr_s1_mc = dFdr_s(x, Af_, x0f_mc_, bf_mc_, cf_mc_);
      double dalphac_dr_s_mc = dFdr_s(x, A_alpha_, x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_);
      // RPA fitting parameters
      double depc_dr_s0_rpa = dFdr_s(x, Ap_, x0p_rpa_, bp_rpa_, cp_rpa_);
      double defc_dr_s1_rpa = dFdr_s(x, Af_, x0f_rpa_, bf_rpa_, cf_rpa_);
      double dalphac_dr_s_rpa = dFdr_s(x, A_alpha_, x0_alpha_rpa_, b_alpha_rpa_,
                                       c_alpha_rpa_);
    
      double fp = two_thirds * (pow((1+zeta),one_third)
             - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      
      // RPA fitting parameters
      double ddelta_e_rpa = defc_dr_s1_rpa - depc_dr_s0_rpa;
      double ddelta_e_mc = defc_dr_s1_mc - depc_dr_s0_mc;
      double ddeltae_rpa_dr_s = f / fpp0 * (1 - zeta4)* dalphac_dr_s_rpa 
                    + f * zeta4 * ddelta_e_mc;
      double ddeltae_rpa_dzeta = alphac_rpa / fpp0 * 
        ( fp * (1.-zeta4) - 4.* f * zeta*zeta2) 
        + delta_e_mc * ( fp*zeta4 + 4.*f*zeta*zeta2);
              
      // Monte Carlo fitting parameters
      // double ddelta_e_mc = defc_dr_s1_mc - depc_dr_s0_mc;
      // double dec_dr_s_mc = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
      //                   + dalphac_dr_s * f / fpp0 * (1 - zeta4);
      // double dec_dzeta_mc = 4.* zeta3 * f * (efc_mc - epc_mc - (alphac_mc/fpp0))
      //        + fp * (zeta4 * (efc_mc - epc_mc) + (1-zeta4)*(alphac_mc/fpp0));
      
      double dec_dr_s = depc_dr_s0_mc + delta_erpa_rszeta/delta_e_rpa * ddelta_e_mc
        + delta_e_mc/delta_e_rpa * ddeltae_rpa_dr_s 
        - delta_erpa_rszeta*delta_e_mc/(delta_e_rpa*delta_e_rpa)* ddelta_e_rpa;
      double dec_dzeta = delta_e_mc / delta_e_rpa * ddeltae_rpa_dzeta;      
 
      od.df_drho_a = ec - (rs/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (rs/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      decrs = dec_dr_s;
      deczeta = dec_dzeta;
    }
}

/////////////////////////////////////////////////////////////////////////////
// VWN5LCFunctional
// Coded by Matt Leininger

#define CLASSNAME VWN5LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public LSDACFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN5LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = LSDACFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN5LCFunctional::VWN5LCFunctional(StateIn& s):
  SavableState(s),
  LSDACFunctional(s)
{
}

VWN5LCFunctional::VWN5LCFunctional()
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
  
}

VWN5LCFunctional::VWN5LCFunctional(const RefKeyVal& keyval):
  LSDACFunctional(keyval)
{
    Ap_  = keyval->doublevalue("Ap", KeyValValuedouble(0.0310907));
    Af_  = keyval->doublevalue("Af", KeyValValuedouble(0.01554535));
    A_alpha_  = keyval->doublevalue("A_alpha", KeyValValuedouble(-1./(6.*M_PI*M_PI)));

    x0p_mc_  = keyval->doublevalue("x0p_mc", KeyValValuedouble(-0.10498));
    bp_mc_   = keyval->doublevalue("bp_mc", KeyValValuedouble(3.72744));
    cp_mc_   = keyval->doublevalue("cp_mc", KeyValValuedouble(12.9352));
    x0f_mc_  = keyval->doublevalue("x0f_mc", KeyValValuedouble(-0.32500));
    bf_mc_   = keyval->doublevalue("bf_mc", KeyValValuedouble(7.06042));
    cf_mc_   = keyval->doublevalue("cf_mc", KeyValValuedouble(18.0578));

    x0_alpha_mc_ = keyval->doublevalue("x0_alpha_mc", KeyValValuedouble(-0.00475840));
    b_alpha_mc_  = keyval->doublevalue("b_alpha_mc", KeyValValuedouble(1.13107));
    c_alpha_mc_  = keyval->doublevalue("c_alpha_mc", KeyValValuedouble(13.0045));

}

VWN5LCFunctional::~VWN5LCFunctional()
{
}

void
VWN5LCFunctional::save_data_state(StateOut& s)
{
  cout << "VWN5LCFunctional: cannot save state" << endl;
  abort();
}

double
VWN5LCFunctional::F(double x, double A, double x0, double b, double c)
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
VWN5LCFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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

#define CLASSNAME XalphaFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
XalphaFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

XalphaFunctional::XalphaFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

XalphaFunctional::XalphaFunctional()
{
  alpha_ = 0.70;
  factor_ = alpha_ * 2.25 * pow(3.0/(4.*M_PI), 1.0/3.0);
}

XalphaFunctional::XalphaFunctional(const RefKeyVal& keyval):
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
  cout << "XalphaFunctional: cannot save state" << endl;
  abort();
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
XalphaFunctional::print(ostream& o) const
{
  o << node0
    << indent << scprintf("XalphaFunctional: alpha = %12.8f", alpha_) << endl;
}

/////////////////////////////////////////////////////////////////////////////
// Becke88XFunctional

#define CLASSNAME Becke88XFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Becke88XFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

Becke88XFunctional::Becke88XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

Becke88XFunctional::Becke88XFunctional()
{
  beta_ = 0.0042;
  beta6_ = 6. * beta_;
  beta26_ = beta6_ * beta_;
  beta2_ = beta_ * beta_;
}

Becke88XFunctional::Becke88XFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  beta_ = keyval->doublevalue("beta", KeyValValuedouble(0.0042));
  beta6_ = 6. * beta_;
  beta26_ = beta6_ * beta_;
  beta2_ = beta_ * beta_;
}

Becke88XFunctional::~Becke88XFunctional()
{
}

void
Becke88XFunctional::save_data_state(StateOut& s)
{
  cout << "Becke88XFunctional: cannot save state" << endl;
  abort();
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
  const double beta26 = beta26_;
  const double beta2 = beta2_;

  // const double beta6=0.0252;
  // const double beta26=0.00010584; 
  // const double beta2=0.00001764;

  // Use simplified formula
  double rho_a_13 = pow(id.a.rho,(1./3.));
  double rho_a_43 = id.a.rho*rho_a_13;
  double xa = sqrt(id.a.gamma)/rho_a_43;
  double xa2 = xa*xa;
  double ga_denom = 1/(1.+beta6*xa*asinh(xa));
  double ga_denom2 = ga_denom*ga_denom;
  double Fa = sqrt(1.+xa2);
  double Ha = 1. - 6.*beta*xa2/Fa;
  double ex = -rho_a_43*beta*xa2*ga_denom;

  if (compute_potential_) {
      od.df_drho_a = 4./3. * beta * rho_a_13 * xa2 * ga_denom2 * Ha;
      od.df_dgamma_aa = -beta * ga_denom / (2.*rho_a_43) * (1. + ga_denom*Ha);
      od.df_drho_b=od.df_drho_a;
      od.df_dgamma_bb=od.df_dgamma_aa;
         }

  if (spin_polarized_) {
      double rho_b_13 = pow(id.b.rho,(1./3.));
      double rho_b_43 = id.b.rho*rho_b_13;
      double xb = sqrt(id.b.gamma)/rho_b_43;
      double xb2 = xb*xb;
      double gb_denom = 1./(1.+beta6*xb*asinh(xb));
      double gb_denom2 = gb_denom*gb_denom;
      double Fb = sqrt(1.+xb2);
      double Hb = 1. - 6.*beta*xb2/Fb;   
      ex += -rho_b_43*beta*xb2*gb_denom;

      if (compute_potential_) {
        od.df_drho_b = 4./3. * beta * rho_b_13 * xb2 * gb_denom2 * Hb;
        od.df_dgamma_bb = -beta / (2.*rho_b_43) * (gb_denom + gb_denom2*Hb);        
      }
    }
  else 
      ex += ex;

  od.energy = ex;
}


/////////////////////////////////////////////////////////////////////////////
// LYPCFunctional
// Coded by Matt Leininger

#define CLASSNAME LYPCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LYPCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

LYPCFunctional::LYPCFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

LYPCFunctional::LYPCFunctional()
{
  a_ = 0.04918;
  b_ = 0.132;
  c_ = 0.2533;
  d_ = 0.349;
}

LYPCFunctional::LYPCFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  a_ = keyval->doublevalue("a", KeyValValuedouble(0.04918));
  b_ = keyval->doublevalue("b", KeyValValuedouble(0.132));
  c_ = keyval->doublevalue("c", KeyValValuedouble(0.2533));
  d_ = keyval->doublevalue("d", KeyValValuedouble(0.349));
}

LYPCFunctional::~LYPCFunctional()
{
}

void
LYPCFunctional::save_data_state(StateOut& s)
{
  cout << "LYPFunctional: cannot save state" << endl;
  abort();
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
       // double ddelta_drho_a = -dens4_3/3.*( c + d/(denom*denom)); 
       double ddelta_drho_a = 1./3* (d*d*dens4_3*dens1_3/(denom*denom)
                                     - delta/dens);
       double domega_drho_a = -1./3.*omega*dens4_3*(11./dens1_3 - c - d/denom);
       //double domega_drho_a = expo/(3.*denom*pow(dens,15./3.) ) *
       //                      ( c + d/denom - 11.*pow(dens,1./3.) );
       double df1_drho_a = -4.*a*dens_ab/(dens*denom)
                         * (1./3.*d*dens4_3/denom + 1/id.a.rho - 1./dens); 
       double  df2_drho_a = -domega_drho_a*a*b*intermediate_3 
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
          double df1_drho_b = -4.*a*dens_ab/(dens*denom)
                            * (1./3.*d*dens4_3/denom + 1./id.b.rho - 1./dens);
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

#define CLASSNAME P86CFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
P86CFunctional::_castdown(const ClassDesc*cd)
 {
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

P86CFunctional::P86CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

P86CFunctional::P86CFunctional()
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

P86CFunctional::P86CFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  a_ = keyval->doublevalue("a", KeyValValuedouble(1.745*0.11));
  C1_ = keyval->doublevalue("C1", KeyValValuedouble(0.001667));
  C2_ = keyval->doublevalue("C2", KeyValValuedouble(0.002568));
  C3_ = keyval->doublevalue("C3", KeyValValuedouble(0.023266));
  C4_ = keyval->doublevalue("C4", KeyValValuedouble(7.389e-6));
  C5_ = keyval->doublevalue("C5", KeyValValuedouble(8.723));
  C6_ = keyval->doublevalue("C6", KeyValValuedouble(0.472));
  C7_ = keyval->doublevalue("C7", KeyValValuedouble(1e4*C4_));
}

P86CFunctional::~P86CFunctional()
{
}

void
P86CFunctional::save_data_state(StateOut& s)
{
  cout << "P86CFunctional: cannot save state" << endl;
  abort();
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
    
  od.energy = fp86;

  if (compute_potential_) {
      double drs_drhoa = -rs/(3.*rho);
      double drs_drhob = drs_drhoa;
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
      od.df_drho_a = dfp86_drhoa;
      od.df_drho_b = dfp86_drhob;

      // gamma part of potential
      // double dPhi_dgamma_aa = Phi/(2.*gamma_total2);
      // double dPhi_dgamma_bb = dPhi_dgamma_aa;
      // double dfp86_dgamma_aa = fp86*(1./gamma_total2 - dPhi_dgamma_aa);
      double prefactor = exp(-Phi)*C_rho/( fzeta*rho43 );
      double dfp86_dgamma_aa = prefactor * (1.-Phi/2.);
      double dfp86_dgamma_bb = dfp86_dgamma_aa;
      od.df_dgamma_aa = dfp86_dgamma_aa;
      od.df_dgamma_bb = dfp86_dgamma_bb;

      // double dPhi_dgamma_ab = 2.*dPhi_dgamma_aa;
      // double dfp86_dgamma_ab = fp86*(2./gamma_total2 - dPhi_dgamma_ab);
      double dfp86_dgamma_ab = prefactor * (2.-Phi);
      od.df_dgamma_ab = dfp86_dgamma_ab;

   }   
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Burke-Ernzerhof (PBE) Correlation Functional
// J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett., 77, 3865, 1996.
// J. P. Perdew, K. Burke, Y. Wang, Phys. Rev. B, 54, 16533, 1996.
// 
// Coded by Matt Leininger

#define CLASSNAME PBECFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PBECFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PBECFunctional::PBECFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PBECFunctional::PBECFunctional()
{
  gamma_ = 0.03109069086965489503494086371273;
  beta_ = 0.06672455060314922;
  local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);
}

PBECFunctional::PBECFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  local_ = keyval->describedclassvalue("local");
  if (local_.null()) local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);
  // in paper
  //gamma_ = keyval->doublevalue("gamma", KeyValValuedouble(0.031091));
  //beta_ = keyval->doublevalue("beta", KeyValValuedouble(0.066725));
  // in PBE.f
  gamma_ = keyval->doublevalue("gamma", KeyValValuedouble(0.03109069086965489503494086371273));
  beta_ = keyval->doublevalue("beta", KeyValValuedouble(0.06672455060314922));
}

PBECFunctional::~PBECFunctional()
{
}

void
PBECFunctional::save_data_state(StateOut& s)
{
  cout << "PBECFunctional: cannot save state" << endl;
  abort();
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

void
PBECFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  double ec_local, dec_local_rs, dec_local_zeta;
  local_->point_lc(id, od, ec_local, dec_local_rs, dec_local_zeta); 

   // Precalculate terms for efficiency
  double rho = id.a.rho + id.b.rho;
  double rs = pow( (3./(4.*M_PI*rho)), (1./3.));
  double zeta = (id.a.rho - id.b.rho)/rho;
  double phi = 0.5*( pow((1.+zeta),(2./3.)) + pow((1.-zeta),(2./3.)) );
  double phi2 = phi*phi;
  double phi3 = phi2*phi;
  double phi4 = phi3*phi;
  double kf = pow( (3.*M_PI*M_PI*rho), (1./3.) );
  double ks = pow( (4.*kf/M_PI), (1./2.) );
  double gamma_aa = id.a.gamma;
  double gamma_bb = id.b.gamma;
  double gamma_ab = id.gamma_ab;
  double gamma_total = sqrt(gamma_aa + gamma_bb + 2.*gamma_ab);
  double t = gamma_total/(2.*ks*phi*rho);
  double t2 = t*t;
  double t3 = t2*t;
  double t4 = t2*t2;
        
  // Compute alpha energy
  double ratio = beta_/gamma_;
  double A = ratio / (exp(-ec_local/(gamma_*phi3)) - 1.);
  double A2 = A*A;
  double q1 = 1.+A*t2;
  double q2 = q1 + A2*t4;
  double X = 1. + ratio*t2*q1/q2;
  double Hpbe = gamma_*phi3*log(X);
  double ec = rho * Hpbe;
  od.energy += ec;

  if (compute_potential_) {
      // d_drhoa part
      double dzeta_drhoa = 1./rho * (1. - zeta);
      double drs_drhoa = -rs/(3.*rho);
      double dec_local_drhoa = dec_local_rs*drs_drhoa + dec_local_zeta*dzeta_drhoa;
      double q3 = pow( (1.+zeta), -1./3.);
      double q3p = pow( (1.-zeta), -1./3.);
      double dphi_drhoa = 1./3. * dzeta_drhoa * (q3 - q3p);
      double dphi3_drhoa = 3.*phi2*dphi_drhoa;
      double dt_drhoa = -t*(dphi_drhoa/phi + 7./6.*1./rho);
      double dt2_drhoa = 2.*t*dt_drhoa;
      double dt4_drhoa = 4.*t3*dt_drhoa;
      double dexp_drhoa = exp(-ec_local/(gamma_*phi3)) *
                          (-dec_local_drhoa/(gamma_*phi3)
                           + 3.*ec_local*dphi_drhoa/(gamma_*phi4));
      double dA_drhoa = -A2/ratio*dexp_drhoa;
      double dA2_drhoa = 2. * A * dA_drhoa;
      double dX_drhoa = ratio*dt2_drhoa*q1/q2 + ratio*t2 *
                        ( (dA_drhoa*t2 + A*dt2_drhoa)/q2 -
                          q1/(q2*q2)*(dA_drhoa*t2+A*dt2_drhoa+dA2_drhoa*t4
                                      + A2*dt4_drhoa) );
      double dHpbe_drhoa = 3.*phi2*gamma_*dphi_drhoa*log(X) + gamma_*phi3/X*dX_drhoa;
      double dfpbe_drhoa = Hpbe + rho*dHpbe_drhoa;
      od.df_drho_a += dfpbe_drhoa;
      
      // d_drhob part
      double dzeta_drhob = 1./rho * (-1. - zeta);
      double drs_drhob = drs_drhoa;
      double dec_local_drhob = dec_local_rs*drs_drhob + dec_local_zeta*dzeta_drhob;
      double dphi_drhob = 1./3. * dzeta_drhob * (q3 - q3p);
      double dphi3_drhob = 3.*phi2*dphi_drhob;
      double dt_drhob = -t*(dphi_drhob/phi + 7./6.*1./rho);
      double dt2_drhob = 2.*t*dt_drhob;
      double dt4_drhob = 4.*t3*dt_drhob;
      double dexp_drhob = exp(-ec_local/(gamma_*phi3)) *
                          (-dec_local_drhob/(gamma_*phi3)
                           + 3.*ec_local*dphi_drhob/(gamma_*phi4));
      double dA_drhob = -A2/ratio*dexp_drhob;
      double dA2_drhob = 2.*A*dA_drhob;
      double dX_drhob = ratio*dt2_drhob*q1/q2 + ratio*t2 *
                        ( (dA_drhob*t2 + A*dt2_drhob)/q2 -
                          q1/(q2*q2)*(dA_drhob*t2+A*dt2_drhob+dA2_drhob*t4
                                      + A2*dt4_drhob) );
      double dHpbe_drhob = 3.*phi2*gamma_*dphi_drhob*log(X) + gamma_*phi3/X*dX_drhob;
      double dfpbe_drhob = Hpbe + rho*dHpbe_drhob;
      od.df_drho_b += dfpbe_drhob;
      
      // d_dgamma_aa part
      double tdt_dgamma_aa = 0.5/( (2.*ks*phi*rho)*(2.*ks*phi*rho) );
      // double dt_dgamma_aa = 0.5*t/(gamma_total*gamma_total);
      // double dt2_dgamma_aa = 2.*t*dt_dgamma_aa;
      // double dt4_dgamma_aa = 4.*t3*dt_dgamma_aa;
      // double dX_tmp =2.*t*ratio*q1/q2 *(1.+ A*t2*(1./q1 - (1.+2.*A*t2)/q2)); 
      // double dX_dgamma_aa = dX_tmp*dt_dgamma_aa;
      // double dX_dgamma_aa = pow( (1./(2.*ks*phi*rho)), 2.)*ratio*q1/q2
      //                      *(1.+ A*t2*(1./q1 - (1.+2.*A*t2)/q2));
      double dX_dgamma_aa = 2.*ratio*tdt_dgamma_aa*q1/q2
                          + ratio*t2*tdt_dgamma_aa*
                            (2.*A/q2 - (2.*A + 4.*A2*t2)*q1/(q2*q2));
      double dHpbe_dgamma_aa = gamma_*phi3/X * dX_dgamma_aa;
      double dfpbe_dgamma_aa = rho*dHpbe_dgamma_aa;
      od.df_dgamma_aa = dfpbe_dgamma_aa;
      
      // d_dgamma_bb part is equal to the aa part for closed and open shell.
      od.df_dgamma_bb = od.df_dgamma_aa;
      
      // d_dgamma_ab part
      double tdt_dgamma_ab = 1./( (2.*ks*phi*rho)*(2.*ks*phi*rho) );
      // double dt_dgamma_ab = t/(gamma_total*gamma_total);
      // double dX_dgamma_ab = dX_tmp*dt_dgamma_ab;
      // double dX_dgamma_ab = 2.*pow( (1./(2.*ks*phi*rho)), 2.)*ratio*q1/q2
      //                      *(1.+ A*t2*(1./q1 - (1.+2.*A*t2)/q2));
      double dX_dgamma_ab = 2.*ratio*tdt_dgamma_ab*q1/q2
                          + ratio*t2*tdt_dgamma_ab*
                            (2.*A/q2 - (2.*A + 4.*A2*t2)*q1/(q2*q2));
      double dHpbe_dgamma_ab = gamma_*phi3/X*dX_dgamma_ab;
      double dfpbe_dgamma_ab = rho*dHpbe_dgamma_ab;
      od.df_dgamma_ab = dfpbe_dgamma_ab;

   }   
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Wang (PW91) Correlation Functional
// J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson,
// and D. J. Singh, Phys. Rev. B, 46, 6671, 1992.
// 
// Coded by Matt Leininger

#define CLASSNAME PW91CFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW91CFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW91CFunctional::PW91CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PW91CFunctional::PW91CFunctional()
{
 local_ = new PW92LCFunctional;
 local_->set_compute_potential(1);
}

PW91CFunctional::PW91CFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  local_ = keyval->describedclassvalue("local");
  if (local_.null()) local_ = new PW92LCFunctional;
  local_->set_compute_potential(1);
}

PW91CFunctional::~PW91CFunctional()
{
}

void
PW91CFunctional::save_data_state(StateOut& s)
{
  cout << "PW91CFunctional: cannot save state" << endl;
  abort();
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
PW91CFunctional::Cxc(double rs)
{
  double a = 23.266;
  double b = 7.389e-3;
  double c = 8.723;
  double d = 0.472;
  double a1 = 2.568;
  double factor = 1e-3;
  double rs2 = rs*rs;

  double res = factor * (a1 + a*rs + b*rs2)/(1. + c*rs + d*rs2 + 10.*b*rs*rs2);
  return res;
  
}

double
PW91CFunctional::dCxc_drho(double rs, double drs_drho, double Cxcrs)
{
  double a = 23.266;
  double b = 7.389e-3;
  double c = 8.723;
  double d = 0.472;
  double a1 = 2.568;
  double factor = 1e-3;
  double rs2 = rs*rs;
  double rs3 = rs2*rs;
  double numer = a1 + a*rs + b*rs2;
  double denom = 1. + c*rs + d*rs2 + 10.*b*rs3;

  double res = factor * drs_drho/denom * ( a + 2.*b*rs - numer/denom * (c + 2.*d*rs + 30.*b*rs2) );
  return res;
}

void
PW91CFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  double ec_local, dec_local_rs, dec_local_zeta;
  local_->point_lc(id, od, ec_local, dec_local_rs, dec_local_zeta); 

   // Precalculate terms for efficiency
  double rho = id.a.rho + id.b.rho;
  double rs = pow( (3./(4.*M_PI*rho)), (1./3.) );
  double zeta = (id.a.rho - id.b.rho)/rho;
  double g = 0.5*( pow((1.+zeta),(2./3.)) + pow((1.-zeta),(2./3.)) );
  double g2 = g*g;
  double g3 = g2*g;
  double g4 = g3*g;
  double kf = pow( (3.*M_PI*M_PI*rho), (1./3.) );
  double ks = pow( (4.*kf/M_PI), (1./2.) );
  double gamma_aa = id.a.gamma;
  double gamma_bb = id.b.gamma;
  double gamma_ab = id.gamma_ab;
  double gamma_total = sqrt(gamma_aa + gamma_bb + 2.*gamma_ab);
  double t = gamma_total/(2.*ks*g*rho);
  double t2 = t*t;
  double t3 = t2*t;
  double t4 = t2*t2;
  double alpha = 0.09;
  double Cc0 = 0.004235;
  double Cx = -0.001667;
  double nu = 16./M_PI * pow( (3.*M_PI*M_PI), (1./3.) );
  double beta = nu*Cc0;
  double beta2 = beta*beta;
  double A_factor = 2.*alpha/beta;
  double ks2 = ks*ks;
  double kf2 = kf*kf;
  
  // Compute alpha energy
  double Aexp = exp(-A_factor*ec_local/(beta*g3));
  double A = A_factor / (Aexp - 1.);
  double A2 = A*A;
  double A3 = A2*A;
  double A4 = A3*A;
  double q1 = t2 + A*t4;
  double q2 = 1. + A*t2 + A2*t4;
  double X = q1/q2;
  double W = 1. + A_factor*X;
  double H0pw91 = beta/A_factor*g3*log(W);
  double Cxcrs = Cxc(rs);
  double Ccrs = Cxcrs - Cx;
  double Z = 100.*g4*ks2/kf2*t2;
  double expH1 = exp(-Z);
  double Y = nu * (Ccrs - Cc0 - 3.*Cx/7.); 
  double H1pw91 = Y*g3*t2*expH1;
  double Hpw91 = H0pw91 + H1pw91;
  double ec = rho * Hpw91;

  od.energy = ec;

  double rs2 = rs*rs;
  double rs3 = rs2*rs;
  double rs4 = rs3*rs;

  if (compute_potential_) {
      // d_drhoa part
      double dzeta_drhoa = 1./rho * (1. - zeta);
      double drs_drhoa = -rs/(3.*rho);
      double dec_local_drhoa = dec_local_rs*drs_drhoa + dec_local_zeta*dzeta_drhoa;
      double q3 = pow( (1.+zeta), -1./3.);
      double q3p = pow( (1.-zeta), -1./3.);
      double dg_drhoa = 1./3. * dzeta_drhoa * (q3 - q3p);
      double dg3_drhoa = 3.*g2*dg_drhoa;
      double dkf_drhoa = -9.*M_PI/4. * drs_drhoa/rs2;
      double dks_drhoa = 2./(M_PI*ks)*dkf_drhoa;
      double dt_drhoa = -t/(2.*g*ks*rho) * (2.*dg_drhoa*ks*rho + 2.*g*dks_drhoa*rho + 2.*g*ks);
      double dt2_drhoa = 2.*t*dt_drhoa;
      double dt4_drhoa = 4.*t3*dt_drhoa;
      double dks2_drhoa = 2.*ks*dks_drhoa;
      double dkf2_drhoa = 2.*kf*dkf_drhoa;
      // double dAexp_drhoa = Aexp * A_factor * ( -dec_local_drhoa/(g3*beta)
      //                       + ec_local/(g3*g3*beta2*beta) * 3.*g2*beta2*dg_drhoa );
      double dAexp_drhoa = -A_factor/beta*Aexp * (dec_local_drhoa/g3 - 3.*ec_local*dg_drhoa/g4);
      double dA_drhoa = -A/(Aexp - 1.) * dAexp_drhoa;
      double dA2_drhoa = 2. * A * dA_drhoa;
      double dX_drhoa = (2.*t*dt_drhoa + dA_drhoa*t4 + 4.*A*t3*dt_drhoa)/q2
         - X/q2 * (dA_drhoa*t2 + 2.*A*t*dt_drhoa + 2.*A*t4*dA_drhoa + 4.*A2*t3*dt_drhoa);
      double dH0pw91_drhoa = 3.*g2*dg_drhoa*beta/A_factor*log(1.+A_factor*X) +
                             g3*beta*dX_drhoa/(1.+A_factor*X);
      double dCcrs_drhoa = dCxc_drho(rs, drs_drhoa, Cxcrs);
      double dZ_drhoa = Z * (4./g*dg_drhoa+2./ks*dks_drhoa-2./kf*dkf_drhoa+2./t*dt_drhoa);
      double dexpH1_drhoa = -dZ_drhoa*expH1;
      double dY_drhoa = nu * dCcrs_drhoa;
      double dH1pw91_drhoa = dY_drhoa*g3*t2*expH1 + Y *
                             (3.*g2*dg_drhoa*t2*expH1 + 2.*g3*t*dt_drhoa*expH1
                              + g3*t2*dexpH1_drhoa);
      double dHpw91_drhoa = dH0pw91_drhoa + dH1pw91_drhoa;
      double dfpw91_drhoa = Hpw91 + rho*dHpw91_drhoa;
      od.df_drho_a = dfpw91_drhoa;
          
      if (spin_polarized_) {
        // d_drhob part
        double dzeta_drhob = 1./rho * (-1. - zeta);
        double drs_drhob = -rs/(3.*rho);
        double dec_local_drhob = dec_local_rs*drs_drhob + dec_local_zeta*dzeta_drhob;
        double q3 = pow( (1.+zeta), -1./3.);
        double q3p = pow( (1.-zeta), -1./3.);
        double dg_drhob = 1./3. * dzeta_drhob * (q3 - q3p);
        double dg3_drhob = 3.*g2*dg_drhob;
        double dkf_drhob = -9.*M_PI/4. * drs_drhob/rs2;
        double dks_drhob = 2./(M_PI*ks)*dkf_drhob;
        double dt_drhob = -t/(2.*g*ks*rho) *
                          (2.*dg_drhob*ks*rho + 2.*g*dks_drhob*rho + 2.*g*ks);
        double dt2_drhob = 2.*t*dt_drhob;
        double dt4_drhob = 4.*t3*dt_drhob;
        double dks2_drhob = 2.*ks*dks_drhob;
        double dkf2_drhob = 2.*kf*dkf_drhob;
        double dAexp_drhob = Aexp * A_factor * ( -dec_local_drhob/(g3*beta)
                                    + ec_local/(g3*g3*beta2*beta) * 3.*g2*beta2*dg_drhob );
        double dA_drhob = -A/(Aexp - 1.) * dAexp_drhob;
        double dA2_drhob = 2. * A * dA_drhob;
        double dX_drhob = (2.*t*dt_drhob + dA_drhob*t4 + 4.*A*t3*dt_drhob)/q2
           - X/q2 * (dA_drhob*t2 + 2.*A*t*dt_drhob + 2.*A*t4*dA_drhob + 4.*A2*t3*dt_drhob);
        double dH0pw91_drhob = 3.*g2*dg_drhob*beta/A_factor*log(1.+A_factor*X) +
                               g3*beta*dX_drhob/(1.+A_factor*X);
        double dCcrs_drhob = dCxc_drho(rs, drs_drhob, Cxcrs);
        double dZ_drhob = Z*(4./g*dg_drhob+2./ks*dks_drhob-2./kf*dkf_drhob+2./t*dt_drhob);
        double dexpH1_drhob = -dZ_drhob*expH1;
        double dY_drhob = nu*dCcrs_drhob;
        double dH1pw91_drhob = dY_drhob*g3*t2*expH1 + Y *
                             (3.*g2*dg_drhob*t2*expH1 + 2.*g3*t*dt_drhob*expH1
                              + g3*t2*dexpH1_drhob);
        double dHpw91_drhob = dH0pw91_drhob + dH1pw91_drhob;
        double dfpw91_drhob = Hpw91 + rho*dHpw91_drhob;
        od.df_drho_b = dfpw91_drhob;
        }      
      else od.df_drho_b = od.df_drho_a;
      
      // d_dgamma_aa part
      double dt_dgamma_aa = 0.5*t/(gamma_total*gamma_total);
      double dt2_dgamma_aa = 2.*t*dt_dgamma_aa;
      double dt4_dgamma_aa = 4.*t3*dt_dgamma_aa;
      double dX_tmp = (1.-A*X)/q2 * (2.*t + 4.*A*t3);
      double dX_dgamma_aa = dX_tmp*dt_dgamma_aa;
      double dH0pw91_dgamma_aa = g3*beta/(1.+A_factor*X) * dX_dgamma_aa;
      // double dH0pw91_dgamma_aa = g3*beta2/(A_factor*X) * dX_dgamma_aa;
      double dZ_dgamma_aa = 2.*Z/t*dt_dgamma_aa;
      double dexpH1_dgamma_aa = -dZ_dgamma_aa* expH1;
      double dH1pw91_dgamma_aa = Y * (g3*2.*t*dt_dgamma_aa*expH1 + g3*t2*dexpH1_dgamma_aa);
      double dHpw91_dgamma_aa = dH0pw91_dgamma_aa + dH1pw91_dgamma_aa;
      double dfpw91_dgamma_aa = rho*dHpw91_dgamma_aa;
      od.df_dgamma_aa = dfpw91_dgamma_aa;
      // d_dgamma_bb part is equal to the aa part for closed and open shell.
      od.df_dgamma_bb = od.df_dgamma_aa;
      
      // d_dgamma_ab part
      double dt_dgamma_ab = t/(gamma_total*gamma_total);
      double dX_dgamma_ab = dX_tmp*dt_dgamma_ab;
      double dH0pw91_dgamma_ab = g3*beta/(1.+A_factor*X) * dX_dgamma_ab;
      double dZ_dgamma_ab = 2.*Z/t*dt_dgamma_ab;
      double dexpH1_dgamma_ab = -dZ_dgamma_ab * expH1;
      double dH1pw91_dgamma_ab = Y *(g3*2.*t*dt_dgamma_ab*expH1 + g3*t2*dexpH1_dgamma_ab);
      double dHpw91_dgamma_ab = dH0pw91_dgamma_ab + dH1pw91_dgamma_ab;
      double dfpw91_dgamma_ab = rho*dHpw91_dgamma_ab;
      od.df_dgamma_ab = dfpw91_dgamma_ab;

   }   
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Burke-Ernzerhof (PBE) Exchange Functional
// J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett., 77, 3865, 1996.
// J. P. Perdew, K. Burke, Y. Wang, Phys. Rev. B, 54, 16533, 1996.
// 
// Coded by Matt Leininger

#define CLASSNAME PBEXFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PBEXFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PBEXFunctional::PBEXFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PBEXFunctional::PBEXFunctional()
{
  mu_ = 0.2195149727645171;
  kappa_ = 0.804;
}

PBEXFunctional::PBEXFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  // in PBE.f
  mu_ = keyval->doublevalue("mu", KeyValValuedouble(0.2195149727645171));
  // in paper
  //mu_ = keyval->doublevalue("mu", KeyValValuedouble(0.21951));
  kappa_ = keyval->doublevalue("kappa", KeyValValuedouble(0.804));
  int revPBEX = keyval->intvalue("revPBEX", KeyValValueint(0));
  if (revPBEX) kappa_ = 1.245;
}

PBEXFunctional::~PBEXFunctional()
{
}

void
PBEXFunctional::save_data_state(StateOut& s)
{
  cout << "PBEXFunctional: cannot save state" << endl;
  abort();
}

int
PBEXFunctional::need_density_gradient()
{
  return 1;
}

void
 PBEXFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double rhoa = 2. * id.a.rho;
  double k_Fa = pow( (3.*M_PI*M_PI*rhoa), (1./3.) );
  double e_xa_unif = -3. * k_Fa / (4.*M_PI);
  double gamma_aa = 2.*sqrt(id.a.gamma);
  double sa = gamma_aa/(2. * k_Fa * rhoa);
  double sa2 = sa*sa;

  double f_xa = 1. + kappa_ - kappa_ / ( 1. + (mu_*sa2/kappa_) );
  double ex = 0.5 * rhoa * e_xa_unif * f_xa;
  if (compute_potential_) {
    double dex_unif_drhoa = e_xa_unif / (3.*id.a.rho);
    double dsa_drhoa = -4./3. * sa / id.a.rho;
    double dFs_drhoa = kappa_ / pow( ( 1. + (mu_*sa2/kappa_)), 2.) *
                       mu_/kappa_ * 2.*sa*dsa_drhoa;
    double dEx_drhoa = 2.*e_xa_unif*f_xa + rhoa*f_xa*dex_unif_drhoa
                     + rhoa*e_xa_unif*dFs_drhoa;
    od.df_drho_a = 0.5 * dEx_drhoa;
    od.df_dgamma_aa = 0.5 * rhoa * e_xa_unif * mu_ * sa2 / id.a.gamma * 
                      pow( (1. + mu_*sa2/kappa_), -2.);
    od.df_drho_b = od.df_drho_a;
    od.df_dgamma_bb = od.df_dgamma_aa;
    od.df_dgamma_ab = 0.;
  }

  if (spin_polarized_) {
    double rhob = 2. * id.b.rho;
    double k_Fb = pow( (3.*M_PI*M_PI*rhob), (1./3.) );
    double e_xb_unif = -3.*k_Fb/(4.*M_PI);
    double gamma_bb = 2.*sqrt(id.b.gamma);
    double sb = gamma_bb/(2.*k_Fb*rhob);
    double sb2 = sb*sb;
    double f_xb = 1. + kappa_ - kappa_ /( 1. + (mu_*sb2/kappa_) );
    ex += 0.5 * rhob * e_xb_unif * f_xb;
    if (compute_potential_) {
      double dex_unif_drhob = e_xb_unif / (3.*id.b.rho);
      double dsb_drhob = -4./3. * sb/id.b.rho;
      double dFs_drhob = kappa_ / pow( ( 1. + (mu_*sb2/kappa_)), 2.) *
                       mu_/kappa_ * 2.*sb*dsb_drhob;
      double dEx_drhob = 2.*e_xb_unif*f_xb + rhob*f_xb*dex_unif_drhob
                     + rhob*e_xb_unif*dFs_drhob;
      od.df_drho_b = 0.5 * dEx_drhob;
      od.df_dgamma_bb = 0.5 * rhob * e_xb_unif * mu_ * sb2 / id.b.gamma *
              pow( (1. + mu_*sb2/kappa_), -2.);
    }
  }
  else ex += ex;

  od.energy = ex;
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Wang (PW91) Exchange Functional
// J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson,
// and D. J. Singh, Phys. Rev. B, 46, 6671, 1992.
// 
// Coded by Matt Leininger

#define CLASSNAME PW91XFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW91XFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW91XFunctional::PW91XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PW91XFunctional::PW91XFunctional()
{
  a1_ = 0.19645;
  a2_ = 7.7956;
  a3_ = 0.2743;
  a4_ = -0.1508;
  a5_ = 0.004;
  b_ = 100.;
}

PW91XFunctional::PW91XFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  a1_ = keyval->doublevalue("a1", KeyValValuedouble(0.19645));
  a2_ = keyval->doublevalue("a2", KeyValValuedouble(7.7956));
  a3_ = keyval->doublevalue("a3", KeyValValuedouble(0.2743));
  a4_ = keyval->doublevalue("a4", KeyValValuedouble(-0.1508));
  a5_ = keyval->doublevalue("a5", KeyValValuedouble(0.004));
  b_  = keyval->doublevalue("b", KeyValValuedouble(100.));
}

PW91XFunctional::~PW91XFunctional()
{
}

void
PW91XFunctional::save_data_state(StateOut& s)
{
  cout << "PW91XFunctional: cannot save state" << endl;
  abort();
}

int
PW91XFunctional::need_density_gradient()
{
  return 1;
}

void
 PW91XFunctional::point(const PointInputData &id,
                         PointOutputData &od)
{
  od.zero();

  double rhoa = 2. * id.a.rho;
  double r_sa = pow( (3./(4.*M_PI*rhoa)), (1./3.) );
  double r_sa2 = r_sa*r_sa;
  double s_factor = pow( ((9.*M_PI)/4.),(1./3.) );
  double k_Fa = s_factor / r_sa;
  double e_xa_unif = -3. * k_Fa / (4.*M_PI);
  double gamma_aa = 2. * sqrt(id.a.gamma);
  double sa = gamma_aa/(2. * k_Fa * rhoa);
  double sa2 = sa*sa;
  double sa3 = sa2*sa;
  double sa4 = sa2*sa2;
  double sinha = asinh(a2_*sa);
  double asinha = a1_*sa*sinha;
  double expa = exp(-b_*sa2);
  double numera = 1. + asinha + (a3_ + a4_*expa)*sa2;
  double denoma = 1. + asinha + a5_*sa4;
  double Fxa = numera/denoma;
  double fpw91xa = rhoa * e_xa_unif * Fxa;
  double ex = 0.5 * fpw91xa;

  if (compute_potential_) {
    double drs_drhoa = -r_sa/(3.*rhoa);
    double dsa_drhoa = sa * (drs_drhoa/r_sa - 1./rhoa);
    double dsinha_drhoa = 1./sqrt(1.+a2_*a2_*sa2) * a2_*dsa_drhoa;
    double dexpa_drhoa = -2.*b_*sa*dsa_drhoa*expa;
    double dFxa_drhoa = 1./denoma * ( a1_*dsa_drhoa*sinha + a1_*sa*dsinha_drhoa
        + 2.*a3_*sa*dsa_drhoa + 2.*a4_*sa*dsa_drhoa*expa + a4_*sa2*dexpa_drhoa)
                      - Fxa/denoma * (a1_*dsa_drhoa*sinha + a1_*sa*dsinha_drhoa
                                     + 4.*a5_*sa3*dsa_drhoa);
    double dfpw91xa_drhoa = -3./(4.*M_PI)*s_factor*( (1./r_sa - rhoa/r_sa2*drs_drhoa)*Fxa
                            + rhoa/r_sa*dFxa_drhoa);
    od.df_drho_a = 0.5 * dfpw91xa_drhoa; 

    double dsa_dgamma_aa = sa/(2.*gamma_aa*gamma_aa);
    double dsinha_dgamma_aa = 1./sqrt(1.+a2_*a2_*sa2) * a2_*dsa_dgamma_aa;
    double dexpa_dgamma_aa = -2.*b_*sa*dsa_dgamma_aa*expa;
    double dFxa_dgamma_aa = 1./denoma * (a1_*dsa_dgamma_aa*sinha + a1_*sa*dsinha_dgamma_aa
         + 2.*a3_*sa*dsa_dgamma_aa + 2.*a4_*sa*dsa_dgamma_aa*expa + a4_*sa2*dexpa_dgamma_aa)
         - Fxa/denoma
         *(a1_*dsa_dgamma_aa*sinha + a1_*sa*dsinha_dgamma_aa + 4.*a5_*sa3*dsa_dgamma_aa);
    double dfpw91xa_dgamma_aa = -3./(4.*M_PI)*s_factor*rhoa/r_sa * dFxa_dgamma_aa;
    od.df_dgamma_aa = 0.5 * dfpw91xa_dgamma_aa;

      od.df_drho_b = od.df_drho_a;
      od.df_dgamma_bb = od.df_dgamma_aa;
      od.df_dgamma_ab = 0.;
    }

#if 0
  if (spin_polarized_) {
      double rhob = 2. * id.b.rho;
      double r_sb = pow( (3./(4.*M_PI*rhob)), (1./3.) );
      double r_sb2 = r_sb*r_sb;
      double k_Fb = s_factor * 1./r_sb;
      double e_xb_unif = -3.*k_Fb/(4.*M_PI);
      double gamma_bb = 2.*sqrt(id.b.gamma);
      double sb = gamma_bb/(2.*k_Fb*rhob);
      double sb2 = sb*sb;
      double sb3 = sb2*sb;
      double sb4 = sb2*sb2;
      double sinhb = asinh(a2_*sb);
      double asinhb = a1_*sb*sinhb;
      double expb = exp(-b_*sb2);
      double numerb = 1. + asinhb + (a3_ + a4_*expb)*sb2;
      double denomb = 1. + asinhb + a5_*sb4;
      double Fxb = numerb/denomb;
      double fpw91xb = rhob * e_xb_unif * Fxb;

      ex += 0.5 * fpw91xb;

      if (compute_potential_) {
        double drs_drhob = -r_sb/(3.*rhob);
        double dsb_drhob = sb * (1./r_sb*drs_drhob - 1./rhob);
        double dsinhb_drhob = 1./sqrt(1.+a2_*a2_*sb2) * a2_*dsb_drhob;
        double dexpb_drhob = -2.*b_*sb*dsb_drhob*expb;
        double dFxb_drhob = 1./denomb * ( a1_*dsb_drhob*sinhb + a1_*sb*dsinhb_drhob
                + 2.*a3_*sb*dsb_drhob + 2.*a4_*sb*dsb_drhob*expb + a4_*sb2*dexpb_drhob)
                      - Fxb/denomb * (a1_*dsb_drhob*sinhb + a1_*sb*dsinhb_drhob
                                     + 4.*a5_*sb3*dsb_drhob);
        double dfpw91xb_drhob = -3./(4.*M_PI)*s_factor*( (1./r_sb - rhob/r_sb2*drs_drhob)*Fxb
                            + rhob/r_sb*dFxb_drhob);
        od.df_drho_b = 0.5 * dfpw91xb_drhob; 

        double dsb_dgamma_bb = sb/(2.*gamma_bb*gamma_bb);
        double dsinhb_dgamma_bb = 1./sqrt(1.+a2_*a2_*sb2) * a2_*dsb_dgamma_bb;
        double dexpb_dgamma_bb = -2.*b_*sb*dsb_dgamma_bb*expb;
        double dFxb_dgamma_bb = 1./denomb * (a1_*dsb_dgamma_bb*sinhb + a1_*sb*dsinhb_dgamma_bb
           + 2.*a3_*sb*dsb_dgamma_bb + 2.*a4_*sb*dsb_dgamma_bb*expb + a4_*sb2*dexpb_dgamma_bb)
           - Fxb/denomb
           *(a1_*dsb_dgamma_bb*sinhb + a1_*sb*dsinhb_dgamma_bb + 4.*a5_*sb3*dsb_dgamma_bb);
        double dfpw91xb_dgamma_bb = -3./(4.*M_PI)*s_factor*rhob/r_sb * dFxb_dgamma_bb;
        od.df_dgamma_bb = 0.5 * dfpw91xb_dgamma_bb;

       }
    }
  else ex += ex;
#endif
  ex += ex;
  od.energy = ex;
}

/////////////////////////////////////////////////////////////////////////////
// Perdew-Wang (PW86) Exchange Functional
// J. P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800, 1986. 
// 
// Coded by Matt Leininger

#define CLASSNAME PW86XFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW86XFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW86XFunctional::PW86XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PW86XFunctional::PW86XFunctional()
{
  m_ = 1./15.;
  a_ = 0.0864/m_;
  b_ = 14.;
  c_ = 0.2;
}

PW86XFunctional::PW86XFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  m_ = keyval->doublevalue("m", KeyValValuedouble(1./15.));
  a_ = keyval->doublevalue("a", KeyValValuedouble(0.0864/m_));
  b_ = keyval->doublevalue("b", KeyValValuedouble(14.));
  c_ = keyval->doublevalue("c", KeyValValuedouble(0.2));
}

PW86XFunctional::~PW86XFunctional()
{
}

void
PW86XFunctional::save_data_state(StateOut& s)
{
  cout << "PW86XFunctional: cannot save state" << endl;
  abort();
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
  double r_sa = pow( (3./(4.*M_PI*rhoa)), (1./3.) );
  double r_sa2 = r_sa*r_sa;
  double s_factor = pow( ((9.*M_PI)/4.),(1./3.) );
  double k_Fa = s_factor * 1./r_sa;
  // double k_Fa = pow( (3.*M_PI*M_PI*rhoa), (1./3.) );
  double rhoa43 = pow(rhoa, (4./3.));
  double rhoa13 = pow(rhoa, (1./3.));
  double Ax = -3./4.*pow( (3./M_PI), (1./3.) );
  double gamma_aa = 2. * sqrt(id.a.gamma);
  double sa = gamma_aa/(2. * k_Fa * rhoa);
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
    double drs_drhoa = -pow( (3./(4.*M_PI)), (1./3.) ) / 3. * pow(id.a.rho, -4./3.);
    double dsa_drhoa = id.a.gamma/( pow(2., 4./3.)*pow( (3.*M_PI*M_PI), 1./3.) ) *
                                    -4./3.* pow(id.a.rho, 1./3.);
    double dFxa_drhoa = m_ * pow(F1a, (m_-1.)) *
                        ( 2.*a_*sa*dsa_drhoa + 4.*b_*sa3*dsa_drhoa + 6.*c_*sa5*dsa_drhoa);
    double dfpw86xa_drhoa = Ax * ( pow(2., 4./3.)*4./3.*pow(id.a.rho,1./3.)*Fxa
                                   + rhoa43 * dFxa_drhoa );
                           
    od.df_drho_a = 0.5 * dfpw86xa_drhoa; 

    double dsa_dgamma_aa = sa/(2.*gamma_aa*gamma_aa);
    double dFxa_dgamma_aa = m_ * pow(F1a, (m_-1.)) *
                    ( 2.*a_*sa*dsa_dgamma_aa + 4.*b_*sa3*dsa_dgamma_aa + 6.*c_*sa5*dsa_dgamma_aa);
    double dfpw86xa_dgamma_aa = Ax * rhoa43 * dFxa_dgamma_aa;
    od.df_dgamma_aa = 0.5 * dfpw86xa_dgamma_aa;

      od.df_drho_b = od.df_drho_a;
      od.df_dgamma_bb = od.df_dgamma_aa;
      od.df_dgamma_ab = 0.;
    }
  
    if (spin_polarized_) {
      double rhob = 2. * id.b.rho;
      double r_sb = pow( (3./(4.*M_PI*rhob)), (1./3.) );
      double r_sb2 = r_sb*r_sb;
      double k_Fb = s_factor * 1./r_sb;
      double rhob43 = pow(rhob, (4./3.));
      double rhob13 = pow(rhob, (1./3.));
      double gamma_bb = 2.*sqrt(id.b.gamma);
      double sb = gamma_bb/(2.*k_Fb*rhob);
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
          double drs_drhob = -r_sb/(3./rhob);
          double dsb_drhob = sb * (1./r_sb*drs_drhob - 1./rhob);
          double dFxb_drhob = m_ * pow(F1b, (m_-1.)) *
                       ( 2.*a_*sb*dsb_drhob + 4.*b_*sb3*dsb_drhob + 6.*c_*sb5*dsb_drhob);
          double dfpw86xb = Ax * ( 4./3. * rhob13 * Fxb + rhob43 * dFxb_drhob );
          od.df_drho_b = 0.5 * dfpw86xb;

          double dsb_dgamma_bb = sb/(2.*gamma_bb*gamma_bb);
          double dFxb_dgamma_bb = m_ * pow(F1b, (m_-1.)) *
                    ( 2.*a_*sb*dsb_dgamma_bb + 4.*b_*sb3*dsb_dgamma_bb + 6.*c_*sb5*dsb_dgamma_bb);
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

#define CLASSNAME G96XFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
G96XFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

G96XFunctional::G96XFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

G96XFunctional::G96XFunctional()
{
  b_ = 1./137.;
}

G96XFunctional::G96XFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  b_ = keyval->doublevalue("b", KeyValValuedouble(1./137.));
}

G96XFunctional::~G96XFunctional()
{
}

void
G96XFunctional::save_data_state(StateOut& s)
{
  cout << "G96XFunctional: cannot save state" << endl;
  abort();
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
  double gg96a = alpha - b_*gamma_aa32/(rhoa*rhoa);
  double fxg96a = rhoa43 * gg96a;
  double ex = fxg96a;

  if (compute_potential_) {

      double dgg96a_drhoa = 2.*b_*gamma_aa32/(rhoa*rhoa*rhoa);
      double dfxg96a_drhoa = 4./3.*rhoa13*gg96a + rhoa43*dgg96a_drhoa;
      od.df_drho_a = dfxg96a_drhoa;
      od.df_drho_b = od.df_drho_a;
      
      double dgg96a_dgamma_aa = -3.*b_ / ( 4.*rhoa*rhoa*sqrt(gamma_aa) );
      double dfxg96a_dgamma_aa = rhoa43 * dgg96a_dgamma_aa;
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
      double gg96b = alpha - b_*gamma_bb32/(rhob*rhob);
      double fxg96b = rhob43 * gg96b;
      ex += fxg96b;
    
      if (compute_potential_) {
          double dgg96b_drhob = 2.*b_*gamma_bb32/(rhob*rhob*rhob);
          double dfxg96b_drhob = 4./3.*rhob13*gg96b + rhob43*dgg96b_drhob;
          od.df_drho_b = dfxg96b_drhob;

          double dgg96b_dgamma_bb = -3.*b_ / ( 4.*rhob*rhob*sqrt(gamma_bb) );
          double dfxg96b_dgamma_bb = rhob43 * dgg96b_dgamma_bb;
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
