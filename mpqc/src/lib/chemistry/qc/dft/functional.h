//
// functional.h --- definition of the dft functional
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

#ifndef _chemistry_qc_dft_functional_h
#define _chemistry_qc_dft_functional_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <math/scmat/vector3.h>
#include <chemistry/qc/wfn/wfn.h>

struct PointInputData {
    enum {X=0, Y=1, Z=2};
    enum {XX=0, YX=1, YY=2, ZX=3, ZY=4, ZZ=5};
    struct SpinData {
        double rho;
        // rho^(1/3)
        double rho_13;

        double del_rho[3];
        // gamma = |del rho|
        double gamma;

        // hessian of rho
        double hes_rho[6];
        // del |del rho|
        double del_gamma[3];
        // del^2 rho
        double lap_rho;
        // (del rho).(del |del rho|)
        double del_rho_del_gamma;
    };
    SpinData a, b;

    const SCVector3 &r;

    PointInputData(const SCVector3& r_): r(r_) {}
};

struct PointOutputData {
    // energy at r
    double energy;

    // derivative of functional wrt density
    double df_drho_a;
    double df_drho_b;

    // derivative of functional wrt density gradient
    double df_dgamma_aa;
    double df_dgamma_bb;
    double df_dgamma_ab;
 
  void zero(){energy=df_drho_a=df_drho_b=df_dgamma_aa=df_dgamma_bb=df_dgamma_ab=0.0;}

};

class DenFunctional: virtual public SavableState {
#   define CLASSNAME DenFunctional
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int spin_polarized_;
    int compute_potential_;
    double a0_;  // for ACM functionals

  public:
    DenFunctional();
    DenFunctional(const RefKeyVal &);
    DenFunctional(StateIn &);
    ~DenFunctional();
    void save_data_state(StateOut &);

    // Set to zero if dens_alpha == dens_beta everywhere.
    // The default is false.
    virtual void set_spin_polarized(int i);
    // Set to nonzero if the potential should be computed.
    // The default is false.
    virtual void set_compute_potential(int i);

    // Must return 1 if the density gradient must also be provided.
    // The default implementation returns 0.
    virtual int need_density_gradient();

    virtual void point(const PointInputData&, PointOutputData&) = 0;

    double a0() const { return a0_; }
};
SavableState_REF_dec(DenFunctional);

/** The NElFunctional computes the number of electrons.
    It is primarily for testing the integrator. */
class NElFunctional: public DenFunctional {
#   define CLASSNAME NElFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    NElFunctional();
    NElFunctional(const RefKeyVal &);
    NElFunctional(StateIn &);
    ~NElFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

/** The SumDenFunctional computes energies and densities
    using the a sum of energy density functions method. */
class SumDenFunctional: public DenFunctional {
#   define CLASSNAME SumDenFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    RefDenFunctional *funcs_;
    double *coefs_;
  public:
    SumDenFunctional();
    SumDenFunctional(const RefKeyVal &);
    SumDenFunctional(StateIn &);
    ~SumDenFunctional();
    void save_data_state(StateOut &);

    void set_spin_polarized(int);
    void set_compute_potential(int);
    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);

    void print(ostream& =cout) const;
};

/** The SlaterXFunctional computes energies and densities
    using the Slater exchange term. */
class SlaterXFunctional: public DenFunctional {
#   define CLASSNAME SlaterXFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    SlaterXFunctional();
    SlaterXFunctional(const RefKeyVal &);
    SlaterXFunctional(StateIn &);
    ~SlaterXFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

/** The VWN5CFunctional computes energies and densities using the
    VWN5 (LSDA) correlation term (from Vosko, Wilk, and Nusair). */
class VWN5CFunctional: public DenFunctional {
#   define CLASSNAME VWN5CFunctional
#   define HAVE_KEYVAL_CTOR
 #   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int vwn_; 
    double A_, x0_, b_, c_;
    double F(double x, double A, double x0, double b, double c);
    double dFdr_s(double x, double A, double x0, double b, double c);
  public:
    VWN5CFunctional();
    VWN5CFunctional(const RefKeyVal &);
    VWN5CFunctional(StateIn &);
    ~VWN5CFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

/** The VWN3CFunctional computes energies and densities using the
    VWN3 (LSDA) correlation term (from Vosko, Wilk, and Nusair). */
class VWN3CFunctional: public DenFunctional {
#   define CLASSNAME VWN3CFunctional
#   define HAVE_KEYVAL_CTOR
 #   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int vwn_; 
    double A_, x0_, b_, c_;
    double F(double x, double A, double x0, double b, double c);
    double dFdr_s(double x, double A, double x0, double b, double c);
  public:
    VWN3CFunctional();
    VWN3CFunctional(const RefKeyVal &);
    VWN3CFunctional(StateIn &);
    ~VWN3CFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

//    The PW91LCFunctional computes energies and densities using the
//    PW91 local (LSDA) correlation term from J. P. Perdew and Y. Wang.
//    Phys. Rev. B, 45, 13244, 1992 
class PW92LCFunctional: public DenFunctional {
#   define CLASSNAME PW92LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double F(double x, double A, double alpha_1, double beta_1, double beta_2, 
             double beta_3, double beta_4, double p);
    double dFdr_s(double x, double A, double alpha_1, double beta_1, double beta_2, 
             double beta_3, double beta_4, double p);
  public:
    PW92LCFunctional();
    PW92LCFunctional(const RefKeyVal &);
    PW92LCFunctional(StateIn &);
    ~PW92LCFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

/** The XalphaFunctional computes energies and densities
    using the Xalpha method. */
class XalphaFunctional: public DenFunctional {
#   define CLASSNAME XalphaFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double alpha_;
    double factor_;
  public:
    XalphaFunctional();
    XalphaFunctional(const RefKeyVal &);
    XalphaFunctional(StateIn &);
    ~XalphaFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);

    void print(ostream& =cout) const;
};

/** The Becke88XFunctional computes energies and densities
    Becke's 1988 exchange functional. */
class Becke88XFunctional: public DenFunctional {
#   define CLASSNAME Becke88XFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double beta_;
    double beta6_;
    double beta26_;
    double beta2_;
  public:
    Becke88XFunctional();
    Becke88XFunctional(const RefKeyVal &);
    Becke88XFunctional(StateIn &);
    ~Becke88XFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

/** The LYPCFunctional computes energies and densities
    using the Lee, Yang, and Parr functional. */
class LYPCFunctional: public DenFunctional {
#   define CLASSNAME LYPCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double a_;
    double b_;
    double c_;
    double d_;
  public:
    LYPCFunctional();
    LYPCFunctional(const RefKeyVal &);
    LYPCFunctional(StateIn &);
    ~LYPCFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

// The Perdue-Burke-Ernzerhof 1996 (PBE) exchange functional
class PBEXFunctional: public DenFunctional {
#   define CLASSNAME PBEXFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double mu_;
    double kappa_;     
  public:
    PBEXFunctional();
    PBEXFunctional(const RefKeyVal &);
    PBEXFunctional(StateIn &);
    ~PBEXFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

#if 0
// The PW91 Correlation Functional computes energies and densities
//    using the Lee, Yang, and Parr functional.
// - MLL - What in the heck does LYP have to do with this? 
// The uses the PW91 correlation functional with local correlation
// functional of PW92.
class PW91CFunctional: public DenFunctional {
#   define CLASSNAME PW91CFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void CORPW91(double RS, double ZET, double G, double T,
                 double UU, double VV, double WW,
                 double EC, double ECRS, double ECZET,
                 double& H, double& DVCUP, double& DVCDN);
    void CORLSD(double RS,double ZET,
                double &EC,
                double &VCUP,double &VCDN,
                double &ECRS,double &ECZET,
                double &ALFC);
    void GCOR(double A,double A1,
              double B1, double B2, double B3, double B4,
              double P, double RS, double &GG, double &GGRS);
  public:
    PW91CFunctional();
    PW91CFunctional(const RefKeyVal &);
    PW91CFunctional(StateIn &);
    ~PW91CFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
