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
};

class DenFunctional: virtual_base public SavableState {
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

/** The LSDAXFunctional computes energies and densities
    using the LSDA exchange term. */
class LSDAXFunctional: public DenFunctional {
#   define CLASSNAME LSDAXFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    LSDAXFunctional();
    LSDAXFunctional(const RefKeyVal &);
    LSDAXFunctional(StateIn &);
    ~LSDAXFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};

#if 1
/** The LSDACFunctional computes energies and densities using the
    LSDA correlation term (from Vosko, Wilk, and Nusair). */
class LSDACFunctional: public DenFunctional {
#   define CLASSNAME LSDACFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double F(double x, double A, double x0, double b, double c);
    double dFdr_s(double x, double A, double x0, double b, double c);
  public:
    LSDACFunctional();
    LSDACFunctional(const RefKeyVal &);
    LSDACFunctional(StateIn &);
    ~LSDACFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
};
#endif

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

/** The Becke88Functional computes energies and densities
    Becke's 1988 exchange functional. */
class Becke88Functional: public DenFunctional {
#   define CLASSNAME Becke88Functional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    Becke88Functional();
    Becke88Functional(const RefKeyVal &);
    Becke88Functional(StateIn &);
    ~Becke88Functional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

/** The LYPFunctional computes energies and densities
    using the Lee, Yang, and Parr functional. */
class LYPFunctional: public DenFunctional {
#   define CLASSNAME LYPFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    LYPFunctional();
    LYPFunctional(const RefKeyVal &);
    LYPFunctional(StateIn &);
    ~LYPFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

#if 0
/** The PW91Functional computes energies and densities
    using the Lee, Yang, and Parr functional. */
class PW91Functional: public DenFunctional {
#   define CLASSNAME PW91Functional
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
    PW91Functional();
    PW91Functional(const RefKeyVal &);
    PW91Functional(StateIn &);
    ~PW91Functional();
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
