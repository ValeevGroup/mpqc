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
  abort();
}

DenFunctional::DenFunctional()
{
  spin_polarized_ = 0;
  compute_potential_ = 0;
}

DenFunctional::DenFunctional(const RefKeyVal& keyval)
{
  spin_polarized_ = 0;
  compute_potential_ = 0;
}

DenFunctional::~DenFunctional()
{
}

void
DenFunctional::save_data_state(StateOut& s)
{
  cout << "DenFunctional: cannot save state" << endl;
  abort();
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
  DenFunctional(s)
  maybe_SavableState(s)
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
  cerr << "NElFunctional: cannot save state" << endl;
  abort();
}

void
NElFunctional::point(const PointInputData &id,
                     PointOutputData &od)
{
  od.energy = id.dens_alpha + id.dens_beta;
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
  DenFunctional(s)
  maybe_SavableState(s)
{
}

SumDenFunctional::SumDenFunctional()
{
}

SumDenFunctional::SumDenFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  int ncoef = keyval->count("coefs");
  int nfunc = keyval->count("funcs");
  if (ncoef != nfunc) {
      cerr << "SumDenFunctional: number of coefs and funcs differ" << endl;
      abort();
    }
  n_ = ncoef;
  coefs_ = new double[n_];
  funcs_ = new RefDenFunctional[n_];
  for (int i=0; i<n_; i++) {
      coefs_[i] = keyval->doublevalue("coefs", i);
      funcs_[i] = keyval->describedclassvalue("funcs", i);
    }
}

SumDenFunctional::~SumDenFunctional()
{
  for (int i=0; i<n_; i++) funcs_[i] = 0; // just in case
  delete[] funcs_;
  delete[] coefs_;
}

void
SumDenFunctional::save_data_state(StateOut& s)
{
  cerr << "SumDenFunctional: cannot save state" << endl;
  abort();
}

int
SumDenFunctional::need_density_gradient()
{
  for (int i=0; i<n_; i++) if (funcs_[i]->need_density_gradient()) return 1;
  return 0;
}

void
SumDenFunctional::set_spin_polarized(int i)
{
  spin_polarized_ = i;
  for (int i=0; i<n_; i++) funcs_[i]->set_spin_polarized(i);
}

void
SumDenFunctional::set_compute_potential(int val)
{
  compute_potential_ = val;
  for (int i=0; i<n_; i++) funcs_[i]->set_compute_potential(val);
}

void
SumDenFunctional::point(const PointInputData &id,
                        PointOutputData &od)
{
  od.energy = 0.0;
  od.alpha_pot = od.beta_pot = 0.0;
  PointOutputData tmpod;
  for (int i=0; i<n_; i++) {
      double e,pa,pb;
      funcs_[i]->point(id, tmpod);
      
      od.energy += coefs_[i] * tmpod.energy;
      if (compute_potential_) {
          od.alpha_pot += coefs_[i] * tmpod.alpha_pot;
          od.beta_pot += coefs_[i] * tmpod.alpha_pot;
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
  DenFunctional(s)
  maybe_SavableState(s)
{
}

XalphaFunctional::XalphaFunctional()
{
  alpha_ = 0.75;
  factor_ = alpha_ * 1.5 * pow(3.0/M_PI, 1.0/3.0);
}

XalphaFunctional::XalphaFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
  alpha_ = 0.75;
  factor_ = alpha_ * 1.5 * pow(3.0/M_PI, 1.0/3.0);
}

XalphaFunctional::~XalphaFunctional()
{
}

void
XalphaFunctional::save_data_state(StateOut& s)
{
  cerr << "XalphaFunctional: cannot save state" << endl;
  abort();
}

void
XalphaFunctional::point(const PointInputData &id,
                        PointOutputData &od)
{
  double density = id.dens_alpha + id.dens_beta;
  od.energy = - factor_ * pow(density, 1.0/3.0) * density * 0.75;

  if (compute_potential_) {
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////
// LSDAXFunctional

#define CLASSNAME LSDAXFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LSDAXFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

LSDAXFunctional::LSDAXFunctional(StateIn& s):
  DenFunctional(s)
  maybe_SavableState(s)
{
}

LSDAXFunctional::LSDAXFunctional()
{
}

LSDAXFunctional::LSDAXFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

LSDAXFunctional::~LSDAXFunctional()
{
}

void
LSDAXFunctional::save_data_state(StateOut& s)
{
  cerr << "LSDAXFunctional: cannot save state" << endl;
  abort();
}

void
LSDAXFunctional::point(const PointInputData &id,
                       PointOutputData &od)
{

  const double mcx2rthird = -0.9305257363491;
  if (!spin_polarized_) {
      od.energy = mcx2rthird * 2.0 * id.dens_alpha * id.dens_alpha13;
    }
  else {
      od.energy = mcx2rthird
                * (id.dens_alpha * id.dens_alpha13
                   +id.dens_beta * id.dens_beta13);
    }

  // this is the same as the above
//   double rho = dens_alpha + dens_beta;
//   double zeta = (dens_alpha - dens_beta)/rho;
//   double f = 1.9236610509315362 * (pow(1.+zeta,4./3.)+pow(1.-zeta,4./3.)-2.);
//   double rs = pow(3./(4.*M_PI*rho), 1./3.);
//   // this was in my notes from the NIST WWW site but looks wrong
//   //double epx = -3.*pow(9./(32.*M_PI*M_PI), 1./3.) / rs;
//   // this seems to be right
//   double epx = -3.*pow(9./(256.*M_PI*M_PI), 1./3.) / rs;
//   double efx = epx * pow(2., 1./3.);
//   double ex = epx + (efx-epx)*f;
//   energy = ex * rho;

  // the same thing again
//   const double cx = 0.7385587663820223;
//   double rho = dens_alpha + dens_beta;
//   double zeta = (dens_alpha - dens_beta)/rho;
//   double f = 1.9236610509315362 * (pow(1.+zeta,4./3.)+pow(1.-zeta,4./3.)-2.);
//   double epx = -cx*pow(rho,1./3.);
//   double efx = epx * pow(2., 1./3.);
//   double ex = epx + (efx-epx)*f;
//   energy = ex * rho;

  if (compute_potential_) {
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////
// LSDACFunctional

#define CLASSNAME LSDACFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LSDACFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

LSDACFunctional::LSDACFunctional(StateIn& s):
  DenFunctional(s)
  maybe_SavableState(s)
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
  cerr << "LSDACFunctional: cannot save state" << endl;
  abort();
}

double
LSDACFunctional::F(double x, double A, double x0, double b, double c)
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

// based on the equations given on a NIST WWW site
void
LSDACFunctional::point(const PointInputData &id,
                       PointOutputData &od)
{
  double rho = id.dens_alpha + id.dens_beta;
  double zeta = (id.dens_alpha - id.dens_beta)/rho;
  double x = pow(3./(4.*M_PI*rho), 1./6.);

  double epc    = F(x, 0.0310907,          -0.10498,    3.72744, 12.9352);
  double efc    = F(x, 0.01554535,         -0.32500,    7.06042, 18.0578);
  // the NIST WWW site gives this
  double alphac = F(x, -1./(6.*M_PI*M_PI), -0.00475840, 1.13107, 13.0045);
  // the MOLPRO WWW docs give this
  //double alphac = F(x, -1./6.*M_PI*M_PI, -0.00475840, 1.13107, 13.0045);

  double f = 1.9236610509315362 * (pow(1.+zeta, 4./3.)+pow(1.-zeta, 4./3.)-2.);
  double fpp0 = 8./9. * 1.9236610509315362;
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double delta_ec_1 = efc - epc;
  double beta = fpp0 * delta_ec_1 / alphac - 1.;
  double delta_ec = alphac * (f/fpp0) * (1. + beta * zeta4);
  double ec = epc + delta_ec;
  od.energy = ec * rho;

  if (compute_potential_) {
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////
// Becke88Functional

#define CLASSNAME Becke88Functional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Becke88Functional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

Becke88Functional::Becke88Functional(StateIn& s):
  DenFunctional(s)
  maybe_SavableState(s)
{
}

Becke88Functional::Becke88Functional()
{
}

Becke88Functional::Becke88Functional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

Becke88Functional::~Becke88Functional()
{
}

void
Becke88Functional::save_data_state(StateOut& s)
{
  cerr << "Becke88Functional: cannot save state" << endl;
  abort();
}

int
Becke88Functional::need_density_gradient()
{
  return 1;
}

// Becke's exchange
// From:  C.W. Murray et al.  Mol. Phys. Vol 78  pp 997-1014 (1993)
// originally coded by Mike Colvin
void
Becke88Functional::point(const PointInputData &id,
                         PointOutputData &od)
{
  double ex;

  // Preset terms from Murray's paper
  double beta=0.0042;

  // Choose between energy functionals based on whether this
  // is a closed shell system or not
  if (!spin_polarized_) {
      // Use simplified formula
      double dens_a4_3=pow(id.dens_alpha,4./3.);
      double x=id.dens_grad_alpha/dens_a4_3;
      ex=-2.*beta*dens_a4_3*x*x/(1.+6.*beta*x*asinh(x));
    }
  else {
      // Use exact formula
      double dens_a4_3=pow(id.dens_alpha,4./3.);
      double x_a=id.dens_grad_alpha/dens_a4_3;
      ex=-beta*dens_a4_3*x_a*x_a/(1.+6.*beta*x_a*asinh(x_a));
      double dens_b4_3=pow(id.dens_beta,4./3.);
      double x_b=id.dens_grad_beta/dens_b4_3;
      ex-=beta*dens_b4_3*x_b*x_b/(1.+6.*beta*x_b*asinh(x_b));
    }

  if (compute_potential_) {
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }

  od.energy = ex;

  //cout << scprintf("B88 = %18.16f", energy) << endl;
}

/////////////////////////////////////////////////////////////////////////////
// LYPFunctional

#define CLASSNAME LYPFunctional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LYPFunctional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

LYPFunctional::LYPFunctional(StateIn& s):
  DenFunctional(s)
  maybe_SavableState(s)
{
}

LYPFunctional::LYPFunctional()
{
}

LYPFunctional::LYPFunctional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

LYPFunctional::~LYPFunctional()
{
}

void
LYPFunctional::save_data_state(StateOut& s)
{
  cerr << "LYPFunctional: cannot save state" << endl;
  abort();
}

int
LYPFunctional::need_density_gradient()
{
  return 1;
}

// Lee-Yang-Parr correlation
// From: Burkhard Miehlich, et al.  Chem Phys. Lett. Vol. 157 pp200-206 (1989)
// originally coded by Mike Colvin
void
LYPFunctional::point(const PointInputData &id,
                     PointOutputData &od)
{
  double ec;

  // Precalculate terms for efficiency
  double dens=id.dens_alpha+id.dens_beta;
  double dens2=dens*dens;
  double grad=id.dens_grad_alpha+id.dens_grad_beta;
  double grad2=grad*grad;
  double dens1_3=pow(dens,-1./3.);

  // Precalculate terms defined in Miehlich's paper
  double a=0.04918;
  double b=0.132;
  double c=0.2533;
  double d=0.347;
  double omega=exp(-c*dens1_3)/(1.+d*dens1_3)*pow(dens,-11./3.);
  double delta=c*dens1_3+d*dens1_3/(1.+d*dens1_3);
  double cf=0.3*pow(3.*M_PI*M_PI,2./3.);

  // Choose between energy functionals based on whether this
  // is a closed shell system or not
  if (!spin_polarized_) {
      // Use simplified formula
      ec= -(a*dens/(1. + d*dens1_3)) +
          a*b*omega*dens2*((3. + 7.*delta)*grad2-
                           72.*cf*pow(dens,8./3.))/72;
    }
  else {
      // Use Miehlich's original formula
      double dens_a2=id.dens_alpha*id.dens_alpha;
      double dens_b2=id.dens_beta*id.dens_beta;
      double grad_a2=id.dens_grad_alpha*id.dens_grad_alpha;
      double grad_b2=id.dens_grad_beta*id.dens_grad_beta;
        
      ec= -4*a*id.dens_alpha*id.dens_beta/
          ((1. + d*dens1_3)*dens) - a*b*omega*
          (-2.*grad2*dens2/3. +
           grad_b2*(2.*dens2/3. - dens_a2) +
           grad_a2*(2.*dens2/3. - dens_b2) +
           id.dens_alpha*id.dens_beta*
           ((47./18. - 7.*delta/18.)*grad2 -
            (5./2. - delta/18.)*(grad_a2 + grad_b2) -
            (-11. + delta)*(grad_a2*id.dens_alpha/dens
                            + grad_b2*id.dens_beta/dens)/9. +
            pow(2.,11./3.)*cf*(pow(id.dens_alpha,8./3.)
                               + pow(id.dens_beta,8./3.))));
    }

  if (compute_potential_) {
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }

  od.energy = ec;
}

/////////////////////////////////////////////////////////////////////////////
// PW91Functional

#define CLASSNAME PW91Functional
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenFunctional
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PW91Functional::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenFunctional::_castdown(cd);
  return do_castdowns(casts,cd);
}

PW91Functional::PW91Functional(StateIn& s):
  DenFunctional(s)
  maybe_SavableState(s)
{
}

PW91Functional::PW91Functional()
{
}

PW91Functional::PW91Functional(const RefKeyVal& keyval):
  DenFunctional(keyval)
{
}

PW91Functional::~PW91Functional()
{
}

void
PW91Functional::save_data_state(StateOut& s)
{
  cerr << "PW91Functional: cannot save state" << endl;
  abort();
}

int
PW91Functional::need_density_gradient()
{
  return 1;
}

// Wang, Perdew, Phys. Rev. B 45, 13244, 1992
// from correspondence from Perdew found on WWW
void
PW91Functional::point(const PointInputData &id,
                      PointOutputData &od)
{
  double EC, VCUP, VCDN, ECRS, ECZET, ALFC;

  double rho = id.dens_alpha + id.dens_beta;
  double zeta = (id.dens_alpha - id.dens_beta)/rho;
  double rs = pow(3./(4.*M_PI*rho), 1./3.);

  CORLSD(rs, zeta, EC, VCUP, VCDN, ECRS, ECZET, ALFC);

  od.energy = EC * rho;

  if (compute_potential_) {
      // not true, really, but abort anyway
      cerr << class_name() << ": cannot compute potential" << endl;
      abort();
    }
}

//  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
//  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
//  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
//     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
//  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
void
PW91Functional::CORLSD(double RS,double ZET,
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
PW91Functional::GCOR(double A,double A1,
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
