//
// am05.cc -- derived from:
// DFunctionals.py Exchange and Correlation Density Functionals
// Translated from python to C++ by Curtis Janssen.
//
// This program is part of the PyQuante quantum chemistry program suite.
//
// Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 
//
// PyQuante version 1.2 and later is covered by the modified BSD
// license. Please see the file LICENSE that is part of this
// distribution.
//

#include <cassert>

#include <util/misc/math.h>

#include <util/misc/scexception.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/dft/am05.h>

using namespace std;
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
norm(double v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}
///////////////////////////////////////////////////////////////////////////
// AM05Functional

static ClassDesc AM05Functional_cd(
  typeid(AM05Functional),"AM05Functional",1,"public DenFunctional",
  0, create<AM05Functional>, create<AM05Functional>);

AM05Functional::AM05Functional(StateIn& s):
  SavableState(s),
  DenFunctional(s)
{
}

AM05Functional::AM05Functional()
{
}

AM05Functional::AM05Functional(const Ref<KeyVal>& keyval):
  DenFunctional(keyval)
{
}

AM05Functional::~AM05Functional()
{
}

void
AM05Functional::save_data_state(StateOut& s)
{
  DenFunctional::save_data_state(s);
}

int
AM05Functional::need_density_gradient()
{
  return 1;
}

int
AM05Functional::need_density_hessian()
{
  return 1;
}

void
AM05Functional::set_spin_polarized(int a)
{
  spin_polarized_ = a;
  if (spin_polarized_) {
      throw FeatureNotImplemented("AM05Functional does not support spin polarized densities", __FILE__, __LINE__, class_desc());
    }
}

void
AM05Functional::point(const PointInputData &id,
                      PointOutputData &od)
{
  // input
  //  n        electron density [bohr**(-3)]
  //  s        scaled abs of gradient of density
  //  u        scaled grad n * grad | grad n |
  //  t        scaled laplacian of density
  //  pot      integer: 1 = calculate potential,
  //             0 = don't calculate potential (then u and t 
  //             are never touched)
  //
  // For exact definitions of s, u, and t, see PRB 33, 8800 (1986).
  // s = |grad(n)|/(2 k_F n)
  // t = (2 k_F)^-2 n^-1 lap(n)
  // u = (2 k_F)^-3 n^-2 grad(n) . grad(|grad(n)|)
  // k_F = (3 pi^2 n)^(1/3)
  double n, s, u, t;
  int pot;

  double rho = id.a.rho + id.b.rho; // Total density 
  double gam = id.a.gamma + id.b.gamma + 2.*id.gamma_ab; // Total gamma
  double fpnt,dfdrho,dfdgamma; // Results
  am05xc(rho,gam,fpnt,dfdrho,dfdgamma);

  od.energy = fpnt;
  od.df_drho_a = dfdrho;
  od.df_drho_b = dfdrho;
  od.df_dgamma_aa =    dfdgamma;
  od.df_dgamma_bb =    dfdgamma;
  od.df_dgamma_ab = 2.*dfdgamma;
}

void
AM05Functional::am05xc(double rho,double gam,
                       double &fxc,double &dfxcdrho,double &dfxcdgamma)
{
  // AM05, both exchange and correlation in one routine
  // AEM June 2006.
  double tol = 1e-16; // to conform with official am05 routine
  fxc = dfxcdrho = dfxcdgamma = 0.0;
  double g = 0.8098;
  double a = 2.804;
  double c = 0.7168;
  if (rho > tol) {
      double snorm2 = (4.*pow(3.*M_PI*M_PI,2./3.)*pow(rho,8./3.));
      double s2 = fabs(gam) / snorm2;
      double s = pow(s2,1./2.);
      // LDAPW exchange and correlation
      double fx0,vxlda;
      xs(0.5*rho,fx0,vxlda); // xs(na) = fx(na) = na*ex(2*na)
      double fxlda = 2.0*fx0;
      double fclda,vclda,vc0b;
      pw(0.5*rho,0.5*rho,fclda,vclda,vc0b);
      assert(vclda == vc0b);
      // Interpolation index
      double X = 1.0/(1.0 + a*s2);
       double w = am05_lambertw(pow(s,3./2.)/sqrt(24.0));
       double zosn;
       if (s < 1.e-14) { // am05_lambertw give back argument if it is < 1.0e-20
           zosn = 1.0; // (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
         }
       else {
          zosn = pow(24.,1./3.)*pow(w,2./3.)/s; // zosn = normalized z/s
        }
       double zfac = s2*pow(zosn*27./32./(M_PI*M_PI),2.0);
       double denom = 1.0 + c*s2*zosn*pow(1.0 + zfac,1./4.); // denom = denominator of Airy LAA refinement function
       double F = (c*s2 + 1.0)/denom; // Airy LAA refinement function
       // Refinement functions
       double Hx = X + (1.0 - X)*F;
       double Hc = X + g*(1.0 - X);
       // Exchange-correlation energy density, Exc = Integrate[fxc]
       fxc = fxlda*Hx + fclda*Hc;
       // Derivatives
       // Interpolation index derivatives: 1/s dX/ds
       double Xsos = -2.0*a*X*X;
       double szsoz = 1.0/(1.0 + w); // szsoz = s*(dz/ds)/z 
       // Airy LAA refinement function derivatives, 1/s dF/ds
       double Fsos = c/(denom*denom)*(2.0 - zosn*
                                      ((1.0 - c*s2)*pow(1.0 + zfac,1./4.) +
                                       (1.0 + c*s2)*(1.0 + 3./2.*zfac)/
                                       pow(1.0 + zfac,3./4.)*szsoz));
       // Refinement function derivatives, 1/s dH{x,c}/ds;
       double Hxsos = (1.0 - X)*Fsos - (F - 1.0)*Xsos;
       double Hcsos = Xsos*(1.0 - g);
       // Exchange-correlation energy density derivatives
       dfxcdrho = vxlda*Hx + vclda*Hc - 4./3.*s2/rho*(Hcsos*fclda +
                                                      Hxsos*fxlda);
       dfxcdgamma = 1./2./snorm2*(Hcsos*fclda + Hxsos*fxlda);
    }
}
    
double
AM05Functional::am05_lambertw(double z)
{
  double result;
  // Used only in am05xc. AEM June, 2006.
  assert(z >= 0.0);
  // If z small, go with the first term of the power expansion, z
  if (z < 1.e-20) {
      result = z;
    }
  else {
      double e = exp(1.0);
      // Inital guess
      if (fabs(z + 1.0/e) > 1.45 ) {
          // Asymptotic expansion at 0 and Inf
          result = log(z);
          result = result - log(result);
        }
       else {
           // Series expansion about -1/e to first order
           result = sqrt(2.0*e*z + 2.0) - 1.0;
         }
      // Find result through iteration
       for (int i=1; i<11; i++) {
           double p = exp(result);
           double t = result*p - z;
           if (result != -1.0) {
               t = t/(p*(result + 1.0) -
                      0.5*(result + 2.0)*t/(result + 1.0));
             }
           else {
              t = 0.0;
            }
           result = result - t;
           if (fabs(t) < (2.48*1.e-14)*(1.0 + fabs(result))) break;
           assert(i != 10);
         }
    }
  return result;
}

void
AM05Functional::xs(double rho,double &ex, double &vx)
{
  // Xalpha X functional
  double tol = 1e-10;
  double Xalpha = 2./3.;
  double fac=-2.25*Xalpha*pow(3./4./M_PI,1./3.);
  if (rho < tol) rho=0;
  double rho3 = pow(rho,1./3.);
  ex = fac*rho*rho3;
  vx = (4./3.)*fac*rho3;
}

void
AM05Functional::pw(double rhoa,double rhob,
                   double &ec,double &vca,double &vcb)
{
  // AEM June 2006.
  double tol = 1e-10;
  double rho = rhoa+rhob;
  ec = vca = vcb = 0;
  if (rho < tol) return;
  double eps;
  cpbe_lsd(rhoa,rhob,eps,vca,vcb);
  ec = rho*eps;
}
    
void
AM05Functional::cpbe_lsd(double rhoa,double rhob,
                         double &eps,double &vca,double &vcb)
{
  // Not quite VWN. AEM: It's usually called PW correlation
  // LSD terms
  // Note that this routine gives out ec, not fc.
  // If you rather have fc, use pw instead
  double rho = rhoa+rhob;
  double Rs = pow(3./(4.*M_PI*rho),1./3.);
  double Zeta = (rhoa-rhob)/rho;
  double thrd = 1./3.;     // thrd*=various multiples of 1/3
  double thrd4 = 4*thrd;
  double ggam=0.5198420997897463295344212145565; // gam= 2^(4/3)-2
  double fzz=8./(9.*ggam); // fzz=f''(0)= 8/(9*gam)
  double rtrs = sqrt(Rs);
  double eu,eurs;
  pbe_gcor(0.0310907,0.21370,7.5957,
           3.5876,1.6382,0.49294,rtrs,
           eu,eurs);
  double ep,eprs;
  pbe_gcor(0.01554535,0.20548,14.1189,
           6.1977,3.3662,0.62517,rtrs,
           ep,eprs);
  double alfm,alfrsm;
  pbe_gcor(0.0168869,0.11125,10.357,
           3.6231,0.88026,0.49671,rtrs,
           alfm,alfrsm);
  double alfc = -alfm;
  double z4 = Zeta*Zeta*Zeta*Zeta;
  double f=(pow(1.+Zeta,thrd4)+pow(1.-Zeta,thrd4)-2.)/ggam;
  eps = eu*(1.-f*z4)+ep*f*z4-alfm*f*(1.-z4)/fzz;

  double ecrs = eurs*(1.-f*z4)+eprs*f*z4-alfrsm*f*(1.-z4)/fzz;
  double fz = thrd4*(pow(1.+Zeta,thrd)-pow(1.-Zeta,thrd))/ggam;
  double eczet = 4.*(Zeta*Zeta*Zeta)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu-(1.-z4)*alfm/fzz);
  double comm = eps -Rs*ecrs/3.-Zeta*eczet;
  vca = comm + eczet;
  vcb = comm - eczet;
}
    
void
AM05Functional::pbe_gcor(double a,double a1,
                         double b1,double b2,double b3,double b4,
                         double rtrs,
                         double &gg,double &ggrs)
{
  //      subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
  // slimmed down version of gcor used in pw91 routines, to interpolate
  // lsd correlation energy, as given by (10) of
  // j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
  // k. burke, may 11, 1996.
  //      implicit real*8 (a-h,o-z)
  double q0 = -2.*a*(1.+a1*rtrs*rtrs);
  double q1 = 2.*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)));
  double q2 = log(1.+1./q1);
  gg = q0*q2;
  double q3 = a*(b1/rtrs+2.*b2+rtrs*(3.*b3+4.*b4*rtrs));
  ggrs = -2.*a*a1*q2-q0*q3/(q1*(1.+q1));
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
