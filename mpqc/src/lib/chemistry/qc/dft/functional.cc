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

#define CLASSNAME PW92LCFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW92LCFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW92LCFunctional::PW92LCFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

PW92LCFunctional::PW92LCFunctional()
{
}

PW92LCFunctional::PW92LCFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
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
  double Q_1 =  2.*A*(beta_1 * x + beta_2 * x2 + beta_3*x*x2 + beta_4 * pow(x2,p));
  double Q_1prime = A * 
           ( beta_1 * 1./x + 2.*beta_2 + 3.*beta_3*x + 2.*(p+1.)*beta_4*pow(x2,p));
  double res = -2.*A*alpha_1*log(1. + 1./Q_1) - Q_0*Q_1prime/(Q_1*Q_1 + Q_1);

  return res;
}

void
PW92LCFunctional::point(const PointInputData &id,
                       PointOutputData &od)
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

  double epc    = F(x, 0.0310907,  0.21370, 7.5957,  3.5876, 1.6382,  0.49294, 1.00);
  double efc    = F(x, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662,  0.62517, 1.00);
  double alphac = F(x, 0.0168869, 0.11125, 10.357,  3.6231, 0.88026, 0.49671, 1.00);
     
  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_ec = -alphac * f / fpp0 * (1. - zeta4) + (efc - epc) * f * zeta4;
  double ec = epc + delta_ec;

  od.energy = ec * rho;

  if (compute_potential_) {
    if (!spin_polarized_) {
      double depc_dr_s0 = 
             dFdr_s(x, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, 1.00);
      double dec_dr_s = depc_dr_s0;
      od.df_drho_a = od.df_drho_b = ec - (x/3.)*dec_dr_s;
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
      od.df_drho_a = ec - (x/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (x/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// VWN5CFunctional

#define CLASSNAME VWN5CFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN5CFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN5CFunctional::VWN5CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

VWN5CFunctional::VWN5CFunctional()
{
}

VWN5CFunctional::VWN5CFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

VWN5CFunctional::~VWN5CFunctional()
{
}

void
VWN5CFunctional::save_data_state(StateOut& s)
{
  cout << "VWN5CFunctional: cannot save state" << endl;
  abort();
}

double
VWN5CFunctional::F(double x, double A, double x0, double b, double c)
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
VWN5CFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
VWN5CFunctional::point(const PointInputData &id,
                       PointOutputData &od)
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

  double epc    = F(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
  double efc    = F(x, 0.01554535,         -0.32500,    7.06042, 18.0578);
  double alphac = F(x, -1./(6.*M_PI*M_PI), -0.00475840, 1.13107, 13.0045);
     
  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double beta = fpp0 * (efc - epc) / alphac - 1.;
  double delta_ec = alphac * f / fpp0 * (1. + beta * zeta4);
  double ec = epc + delta_ec;

  od.energy = ec * rho;

  if (compute_potential_) {
    if (!spin_polarized_) {
      double depc_dr_s0 = dFdr_s(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
      double dec_dr_s = depc_dr_s0;
      od.df_drho_a = od.df_drho_b = ec - (x/3.)*dec_dr_s;
    }
    else {
      double zeta3 = zeta2*zeta;
      double depc_dr_s0 = dFdr_s(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
      double defc_dr_s1 = dFdr_s(x, 0.01554535,         -0.32500,    7.06042, 18.0578);
      double dalphac_dr_s = dFdr_s(x, -1./(6.*M_PI*M_PI), -0.00475840, 1.13107, 13.0045);
      double dec_dr_s = depc_dr_s0*(1 - f*zeta4) + defc_dr_s1 * f * zeta4
                        + dalphac_dr_s * f / fpp0 * (1 - zeta4);
      double fp = two_thirds * (pow((1+zeta),one_third) 
            - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      double dec_dzeta = 4.* zeta3 * f * (efc - epc - (alphac/fpp0))
              + fp * (zeta4 * (efc - epc) + (1-zeta4)*(alphac/fpp0));
      od.df_drho_a = ec - (x/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (x/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
      } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// VWN3CFunctional

#define CLASSNAME VWN3CFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VWN3CFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

VWN3CFunctional::VWN3CFunctional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

VWN3CFunctional::VWN3CFunctional()
{
}

VWN3CFunctional::VWN3CFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

VWN3CFunctional::~VWN3CFunctional()
{
}

void
VWN3CFunctional::save_data_state(StateOut& s)
{
  cout << "VWN3CFunctional: cannot save state" << endl;
  abort();
}

double
VWN3CFunctional::F(double x, double A, double x0, double b, double c)
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
VWN3CFunctional::dFdr_s(double x, double A, double x0, double b, double c)
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
VWN3CFunctional::point(const PointInputData &id,
                       PointOutputData &od)
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
 
  // Monte Carlo fitting parameters 
  double epc_mc    = F(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
  double efc_mc    = F(x, 0.01554535,         -0.32500,    7.06042, 18.0578);
  double alphac_mc = F(x, -1./(6.*M_PI*M_PI), -0.00475840, 1.13107, 13.0045);
  // RPA fitting parameters
  double epc_rpa    = F(x, 0.0310907,          -0.409286,  13.0720,  42.7198);
  double efc_rpa    = F(x, 0.01554535,         -0.743294,  20.1231, 101.578);
  double alphac_rpa = F(x, -1./(6.*M_PI*M_PI), -0.228344,   1.06835, 11.4813);

  double f = 9./8.*fpp0*(pow(1.+zeta, four_thirds)+pow(1.-zeta, four_thirds)-2.);
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_e_rpa = efc_rpa - epc_rpa;
  double delta_e_mc  = efc_mc - epc_mc;
  double beta = fpp0 * delta_e_rpa / alphac_rpa - 1.;
  double delta_erpa_rszeta = alphac_rpa * f / fpp0 * (1. + beta * zeta4);
  double delta_ec = delta_e_mc/delta_e_rpa * delta_erpa_rszeta;
  double ec = epc_rpa + delta_ec;

  od.energy = ec * rho;

  if (compute_potential_) {
      double zeta3 = zeta2*zeta;
      // Monte Carlo fitting parameters
      double depc_dr_s0_mc = dFdr_s(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
      double defc_dr_s1_mc = dFdr_s(x, 0.01554535,         -0.32500,    7.06042, 18.0578);
      double dalphac_dr_s_mc = dFdr_s(x, -1./(6.*M_PI*M_PI), -0.00475840, 1.13107, 13.0045);
      // RPA fitting parameters
      double depc_dr_s0_rpa = dFdr_s(x, 0.0310907,          -0.409286,  13.0720,  42.7198);
      double defc_dr_s1_rpa = dFdr_s(x, 0.01554535,         -0.743294,  20.1231, 101.578);
      double dalphac_dr_s_rpa = dFdr_s(x, -1./(6.*M_PI*M_PI), -0.228344, 1.06835, 11.4813 );
    
      double fp = two_thirds * (pow((1+zeta),one_third)
             - pow((1-zeta),one_third))/(pow(2.,one_third)-1);
      
      // RPA fitting parameters
      double ddelta_e_rpa = defc_dr_s1_rpa - depc_dr_s0_rpa;
      double ddeltae_rpa_dr_s = f / fpp0 * (1 - zeta4)* dalphac_dr_s_rpa 
                    + f * zeta4 * ddelta_e_rpa;
      double ddeltae_rpa_dzeta = dalphac_dr_s_rpa / fpp0 * 
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
 
      od.df_drho_a = ec - (x/3.)*dec_dr_s - (zeta-1)*dec_dzeta;
      od.df_drho_b = ec - (x/3.)*dec_dr_s - (zeta+1)*dec_dzeta;
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
  double xa = id.a.gamma/rho_a_43;
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
      double xb = id.b.gamma/rho_b_43;
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
  double grad_a2=id.a.gamma*id.a.gamma;
  double grad_b2=id.b.gamma*id.b.gamma;
  double grad_ab=id.a.del_rho[0]*id.b.del_rho[0] + id.a.del_rho[1]*id.b.del_rho[1] 
                   + id.a.del_rho[2]*id.b.del_rho[2];

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
       double ddelta_drho_a = -dens4_3/3.*( c + d/(denom*denom)); 
       double domega_drho_a = expo/(3.*denom*pow(dens,15./3.) ) *
                             ( c + d/denom - 11.*pow(dens,1./3.) );
       double df1_drho_a = -4.*a*id.b.rho/(dens*denom)
                * (1.-(id.a.rho*(1.+2./3.*d*dens1_3))/(dens*denom)); 
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
          double domega_drho_b = domega_drho_b;
          double df1_drho_b = -4.*id.a.rho/(dens*denom)
                   * (1.-(id.b.rho*(1.+2./3.*d*dens1_3))/(dens*denom));
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

#if 0
/////////////////////////////////////////////////////////////////////////////
// PW91CFunctional

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
}

PW91CFunctional::PW91CFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
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

// Wang, Perdew, Phys. Rev. B 45, 13244, 1992
// from correspondence from Perdew found on WWW
void
PW91CFunctional::point(const PointInputData &id,
                      PointOutputData &od)
{
  od.zero();
  double EC, VCUP, VCDN, ECRS, ECZET, ALFC;
  double H, DVCUP, DVCDN;

  double rho = id.dens_alpha + id.dens_beta;
  double zeta = (id.dens_alpha - id.dens_beta)/rho;
  double rs = pow(3./(4.*M_PI*rho), 1./3.);

  CORLSD(rs, zeta, EC, VCUP, VCDN, ECRS, ECZET, ALFC);

  od.energy = EC * rho;

  double thrd = 1./3.;
  double thrd2 = 2*thrd;
  double alpha = pow(9*M_PI/4., thrd);
  double pi32 = 3.*M_PI*M_PI;
  double g = 0.5*(pow(1.+zeta, thrd2) + pow(1.-zeta, thrd2));
  double fk = pow(pi32*rho, thrd);
  double sk = sqrt(4.*fk/M_PI);
  double twoksg = 2.*sk*g;
  double agr = sqrt(id.dens_grad_alpha*id.dens_grad_alpha +
                    id.dens_grad_beta*id.dens_grad_beta);
  double t = agr/(twoksg*rho);
  double uu,vv,ww; // only for potential

  CORPW91(rs, zeta, g, t, uu, vv, ww, EC, ECRS, ECZET, H, DVCUP, DVCDN);

  od.energy += H * rho;
  
  if (compute_potential_) {
      // not true, really, but abort anyway
      cout << class_name() << ": cannot compute potential" << endl;
      abort();
    }
}

//  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
//  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
//  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
//     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
//  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
void
PW91CFunctional::CORLSD(double RS,double ZET,
                       double &EC,
                       double &VCUP,double &VCDN,
                       double &ECRS,double &ECZET,
                       double &ALFC)
{
  const double GAM = 0.5198421;
  const double FZZ = 1.709921;
  const double THRD = 0.333333333333;
  const double THRD4 = 1.333333333333;
  double F = (pow(1.+ZET,THRD4)+pow(1.-ZET,THRD4)-2.)/GAM;
  double EU, EURS, EP, EPRS, ALFM, ALFRSM;
  GCOR(0.0310907,0.21370,7.5957,3.5876,1.6382,
       0.49294,1.00,RS,EU,EURS);
  GCOR(0.01554535,0.20548,14.1189,6.1977,3.3662,
       0.62517,1.00,RS,EP,EPRS);
  GCOR(0.0168869,0.11125,10.357,3.6231,0.88026,
       0.49671,1.00,RS,ALFM,ALFRSM);
//  ALFM IS MINUS THE SPIN STIFFNESS ALFC
  ALFC = -ALFM;
  double Z4 = ZET*ZET*ZET*ZET;
  EC = EU*(1.-F*Z4)+EP*F*Z4-ALFM*F*(1.-Z4)/FZZ;
//  ENERGY DONE. NOW THE POTENTIAL:
  if (compute_potential_) {
      ECRS = EURS*(1.-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.-Z4)/FZZ;
      double FZ = THRD4*(pow(1.+ZET,THRD)-pow(1.-ZET,THRD))/GAM;
      ECZET = 4.*(ZET*ZET*ZET)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
                                                      -(1.-Z4)*ALFM/FZZ);
      double COMM = EC -RS*ECRS/3.-ZET*ECZET;
      VCUP = COMM + ECZET;
      VCDN = COMM - ECZET;
    }
}

//  CALLED BY SUBROUTINE CORLSD
void
PW91CFunctional::GCOR(double A,double A1,
                     double B1, double B2, double B3, double B4,
                     double P, double RS, double &GG, double &GGRS)
{
  double P1 = P + 1.;
  double Q0 = -2.*A*(1.+A1*RS);
  double RS12 = sqrt(RS);
  double RS32 = RS12*RS12*RS12;
  double RSP = pow(RS,P);
  double Q1 = 2.*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP);
  double Q2 = log(1.+1./Q1);
  GG = Q0*Q2;
  double Q3 = A*(B1/RS12+2.*B2+3.*B3*RS12+2.*B4*P1*RSP);
  GGRS = -2.*A*A1*Q2-Q0*Q3/(Q1*Q1+Q1);
}

//  GGA91 CORRELATION
//  INPUT RS: SEITZ RADIUS
//  INPUT ZET: RELATIVE SPIN POLARIZATION
//  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
//  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
//  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
//  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
//  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
//  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
void
PW91CFunctional::CORPW91(double RS, double ZET, double G, double T,
                        double UU, double VV, double WW,
                        double EC, double ECRS, double ECZET,
                        double& H, double& DVCUP, double& DVCDN)
{
  H=0;
  if (fabs(EC) < 1e-14)
    return;

  double XNU=15.75592, CC0=0.004235, CX=-0.001667212;
  double ALF=0.09;
  double C1=0.002568, C2=0.023266, C3=7.389e-6, C4=8.723;
  double C5=0.472, C6=7.389e-2, A4=100.;
  double THRDM=1./3., THRD2=2./3.;

  double BET = XNU*CC0;
  double DELT = 2.*ALF/BET;
  double G3 = G*G*G;
  double G4 = G3*G;
  double PON = -DELT*EC/(G3*BET);
  double B = DELT/(exp(PON)-1.);
  double B2 = B*B;
  double T2 = T*T;
  double T4 = T2*T2;
  double T6 = T4*T2;
  double RS2 = RS*RS;
  double RS3 = RS2*RS;
  double Q4 = 1.+B*T2;
  double Q5 = 1.+B*T2+B2*T4;
  double Q6 = C1+C2*RS+C3*RS2;
  double Q7 = 1.+C4*RS+C5*RS2+C6*RS3;
  double CC = -CX + Q6/Q7;
  double R0 = 0.663436444*RS;
  double R1 = A4*R0*G4;
  double COEFF = CC-CC0-3.*CX/7.;
  double R2 = XNU*COEFF*G3;
  double R3 = exp(-R1*T2);
  double H0 = G3*(BET/DELT)*log(1.+DELT*Q4*T2/Q5);
  double H1 = R3*R2*T2;
  H = H0+H1;
}
#endif

/////////////////////////////////////////////////////////////////////////////
// Perdew-Burke-Ernzerhof (PBE) Exchange Functional
// J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett., 77, 3865, 1996.
// J. P. Perdew, K. Burke, Y. Wang, Phys. Rev. B, 54, 16533, 1996.
// 
// Codes by Matt Leininger

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
  mu_ = 0.21951;
  kappa_ = 0.804;
}

PBEXFunctional::PBEXFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  mu_ = keyval->doublevalue("mu", KeyValValuedouble(0.21951));
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
  double gamma_aa = 2. * id.a.gamma;
  double sa = id.a.gamma/(2. * k_Fa * rhoa);
  double sa2 = sa*sa;

  double f_xa = 1. + kappa_ - kappa_ / ( 1. + (mu_*sa2/kappa_) );
  double ex = 0.5 * rhoa * e_xa_unif * f_xa;
  if (compute_potential_) {
    od.df_drho_a = 0.5 * 4./3. * e_xa_unif * 
         ( f_xa - 2. * mu_ * sa2 * pow((1. + (mu_*sa2)/kappa_),-2.));
    od.df_dgamma_aa = -0.5 * rhoa * e_xa_unif * mu_ * sa2 / (id.a.gamma*id.a.gamma) * 
                      pow( (1. + mu_*sa2/kappa_), -2.);
    od.df_drho_b = od.df_drho_a;
    od.df_dgamma_bb = od.df_dgamma_aa;
    od.df_dgamma_ab = 0.;
  }

  if (spin_polarized_) {
    double rhob = 2. * id.b.rho;
    double k_Fb = pow( (3.*M_PI*M_PI*rhob), (1./3.) );
    double e_xb_unif = -3.*k_Fb/(4.*M_PI);
    double gamma_bb = 2.*id.b.gamma;
    double sb = id.b.gamma/(2.*k_Fb*rhob);
    double sb2 = sb*sb;
    double f_xb = 1. + kappa_ - kappa_ /( 1. + (mu_*sb2/kappa_) );
    ex += 0.5 * rhob * e_xb_unif * f_xb;
    if (compute_potential_) {
      od.df_drho_b = 0.5 * 4./3. * e_xb_unif * 
              ( f_xb - 2.*mu_*sb2 * pow((1. + (mu_*sb2)/kappa_),-2.));
      od.df_dgamma_bb = -0.5 * rhob * e_xb_unif * mu_ * sb2 / (id.b.gamma*id.b.gamma) *
              pow( (1. + mu_*sb2/kappa_), -2.);
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
