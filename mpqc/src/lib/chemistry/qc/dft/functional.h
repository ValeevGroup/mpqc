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

#define MIN_DENSITY 10.*DBL_EPSILON

struct PointInputData {
    enum {X=0, Y=1, Z=2};
    enum {XX=0, YX=1, YY=2, ZX=3, ZY=4, ZZ=5};
    struct SpinData {
        double rho;
        // rho^(1/3)
        double rho_13;

        double del_rho[3];
        // gamma = (del rho).(del rho)
        double gamma;

        // hessian of rho
        double hes_rho[6];
        // del^2 rho
        double lap_rho;
    };
    SpinData a, b;

    // gamma_ab = (del rho_a).(del rho_b)
    double gamma_ab;

    const SCVector3 &r;

    // fill in derived quantities
    void compute_derived(int spin_polarized);

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

    void do_fd_point(PointInputData&id,double&in,double&out,
                     double lower_bound, double upper_bound);
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
    void gradient(const PointInputData&, PointOutputData&,
                  double *gradient, int acenter,
                  const RefGaussianBasisSet &basis,
                  const double *dmat_a, const double *dmat_b,
                  const double *bs_values, const double *bsg_values,
                  const double *bsh_values);

    double a0() const { return a0_; }

    void fd_point(const PointInputData&, PointOutputData&);
    void test(const PointInputData &);
    void test();
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

/** The StdDenFunctional provides a shorthand method to construct
    the standard density functionals.  */
class StdDenFunctional: public SumDenFunctional {
#   define CLASSNAME StdDenFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    char *name_;
    void init_arrays(int n);
  public:
    StdDenFunctional();
    /** This KeyVal CTOR does not call the parent's KeyVal CTOR.  The
        "name" keyword is read from the input and is used to initialize the
        functional.  The name must be one of HFK, XALPHA, HFS, HFB, HFG96,
        BLYP, B3LYP, PW91, or PBE. */
    StdDenFunctional(const RefKeyVal &);
    StdDenFunctional(StateIn &);
    ~StdDenFunctional();
    void save_data_state(StateOut &);

    void print(ostream& =cout) const;
};

// The LSDACFunctional computes energies and densities
//    using the designated local correlation functional.
class LSDACFunctional: public DenFunctional {
#   define CLASSNAME LSDACFunctional
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
  public:
    LSDACFunctional();
    LSDACFunctional(const RefKeyVal &);
    LSDACFunctional(StateIn &);
    ~LSDACFunctional();
    void save_data_state(StateOut &);

    void point(const PointInputData&, PointOutputData&);
    virtual 
      void point_lc(const PointInputData&, PointOutputData&, 
                    double &ec_local, double &decrs, double &deczeta) = 0;

};
SavableState_REF_dec(LSDACFunctional);

// The Perdew-Burke-Ernzerhof (PBE) Correlation Functional
// computes energies and densities using the designated
// local correlation functional.
class PBECFunctional: public DenFunctional {
#   define CLASSNAME PBECFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>  
  protected:
    RefLSDACFunctional local_;
    double gamma_;
    double beta_;
  public:
    PBECFunctional();
    PBECFunctional(const RefKeyVal &);
    PBECFunctional(StateIn &);
    ~PBECFunctional();
    void save_data_state(StateOut &);
    int need_density_gradient();
    void point(const PointInputData&, PointOutputData&);
    void set_spin_polarized(int);
  
};

// The Perdew-Wang 1991 Correlation Functional computes energies and densities
//    using the designated local correlation functional.
class PW91CFunctional: public DenFunctional {
#   define CLASSNAME PW91CFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>  
  protected:
    RefLSDACFunctional local_;
  public:
    PW91CFunctional();
    PW91CFunctional(const RefKeyVal &);
    PW91CFunctional(StateIn &);
    ~PW91CFunctional();
    void save_data_state(StateOut &);
    int need_density_gradient();
    double Cxc(double rs);
    double dCxc_drho(double rs, double drs_drho);

    void point(const PointInputData&, PointOutputData&);
    void set_spin_polarized(int);
  
};

// The Perdew 1986 (P86) Correlation Functional computes energies and densities
//    using the designated local correlation functional.
class P86CFunctional: public DenFunctional {
#   define CLASSNAME P86CFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>  
  protected:
    double a_;
    double C1_;
    double C2_;
    double C3_;
    double C4_;
    double C5_;
    double C6_;
    double C7_;
  public:
    P86CFunctional();
    P86CFunctional(const RefKeyVal &);
    P86CFunctional(StateIn &);
    ~P86CFunctional();
    void save_data_state(StateOut &);
    int need_density_gradient();
    void point(const PointInputData&, PointOutputData&);
  
};

// The SlaterXFunctional computes energies and densities
// using the Slater exchange term. 
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

/** The VWNLCFunctional class from which the various VWN
    (Vosko, Wilk and Nusair) local correlation functionals
    (1, 2, 3, 4, 5) are derived. */
class VWNLCFunctional: public LSDACFunctional {
#   define CLASSNAME VWNLCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double Ap_, Af_, A_alpha_;
    double x0p_mc_, bp_mc_, cp_mc_, x0f_mc_, bf_mc_, cf_mc_;
    double x0p_rpa_, bp_rpa_, cp_rpa_, x0f_rpa_, bf_rpa_, cf_rpa_;
    double x0_alpha_mc_, b_alpha_mc_, c_alpha_mc_, x0_alpha_rpa_, b_alpha_rpa_, c_alpha_rpa_;
    double F(double x, double A, double x0, double b, double c);
    double dFdr_s(double x, double A, double x0, double b, double c);
  public:
    VWNLCFunctional();
    VWNLCFunctional(const RefKeyVal &);
    VWNLCFunctional(StateIn &);
    ~VWNLCFunctional();
    void save_data_state(StateOut &);

    virtual
      void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};
/** The VWNTestLCFunctional computes energies and densities using the
    VWNTest local correlation term (from Vosko, Wilk, and Nusair). */
class VWNTestLCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWNTestLCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int monte_carlo_prefactor_;
  public:
    VWNTestLCFunctional();
    VWNTestLCFunctional(const RefKeyVal &);
    VWNTestLCFunctional(StateIn &);
    ~VWNTestLCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};
    
/** The VWN1LCFunctional computes energies and densities using the
    VWN1 local correlation term (from Vosko, Wilk, and Nusair). */
class VWN1LCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWN1LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double Ap_, x0p_, bp_, cp_, Af_, x0f_, bf_, cf_;
  public:
    VWN1LCFunctional();
    VWN1LCFunctional(const RefKeyVal &);
    VWN1LCFunctional(StateIn &);
    ~VWN1LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};

/** The VWN2LCFunctional computes energies and densities using the
    VWN2 local correlation term (from Vosko, Wilk, and Nusair). */
class VWN2LCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWN2LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    VWN2LCFunctional();
    VWN2LCFunctional(const RefKeyVal &);
    VWN2LCFunctional(StateIn &);
    ~VWN2LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};


/** The VWN3LCFunctional computes energies and densities using the
    VWN3 local correlation term (from Vosko, Wilk, and Nusair). */
class VWN3LCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWN3LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int monte_carlo_prefactor_;
   public:
    VWN3LCFunctional();
    VWN3LCFunctional(const RefKeyVal &);
    VWN3LCFunctional(StateIn &);
    ~VWN3LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};

/** The VWN4LCFunctional computes energies and densities using the
    VWN4 local correlation term (from Vosko, Wilk, and Nusair). */
class VWN4LCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWN4LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
   public:
    VWN4LCFunctional();
    VWN4LCFunctional(const RefKeyVal &);
    VWN4LCFunctional(StateIn &);
    ~VWN4LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};

/** The VWN5LCFunctional computes energies and densities using the
    VWN5 local correlation term (from Vosko, Wilk, and Nusair). */
class VWN5LCFunctional: public VWNLCFunctional {
#   define CLASSNAME VWN5LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
  public:
    VWN5LCFunctional();
    VWN5LCFunctional(const RefKeyVal &);
    VWN5LCFunctional(StateIn &);
    ~VWN5LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};

//    The PW92LCFunctional computes energies and densities using the
//    PW92 local (LSDA) correlation term from J. P. Perdew and Y. Wang.
//    Phys. Rev. B, 45, 13244, 1992 
//    This local correlation functional is used in PW91 and PBE.
class PW92LCFunctional: public LSDACFunctional {
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

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
};

//    The PZ81LCFunctional computes energies and densities using the
//    PZ81 local (LSDA) correlation term from
//    J. P. Perdew and A. Zunger, Phys. Rev. B, 23, 5048, 1981.
//    This local correlation functional is used in P86.
class PZ81LCFunctional: public LSDACFunctional {
#   define CLASSNAME PZ81LCFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double Fec_rsgt1(double rs, double beta_1, double beta_2, double gamma);
    double dFec_rsgt1_drho(double rs, double beta_1, double beta_2, double gamma,
                           double &dec_drs);
    double Fec_rslt1(double rs, double A, double B, double C, double D);
    double dFec_rslt1_drho(double rs, double A, double B, double C, double D,
                           double &dec_drs);
  public:
    PZ81LCFunctional();
    PZ81LCFunctional(const RefKeyVal &);
    PZ81LCFunctional(StateIn &);
    ~PZ81LCFunctional();
    void save_data_state(StateOut &);

    void point_lc(const PointInputData&, PointOutputData&, double &, double &, double &);
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

// The Perdew-Wang 1986 (PW86) Exchange functional
class PW86XFunctional: public DenFunctional {
#   define CLASSNAME PW86XFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double a_;
    double b_;
    double c_;
    double m_;
  public:
    PW86XFunctional();
    PW86XFunctional(const RefKeyVal &);
    PW86XFunctional(StateIn &);
    ~PW86XFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

// The Perdew-Wang 1991 (PW91) Exchange functional
class PW91XFunctional: public DenFunctional {
#   define CLASSNAME PW91XFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double a1_;
    double a2_;
    double a3_;
    double a4_;
    double a5_;
    double b_;
  public:
    PW91XFunctional();
    PW91XFunctional(const RefKeyVal &);
    PW91XFunctional(StateIn &);
    ~PW91XFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

// The Perdew-Burke-Ernzerhof 1996 (PBE) exchange functional
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

// The Gill 1996 (G96) exchange functional
class G96XFunctional: public DenFunctional {
#   define CLASSNAME G96XFunctional
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double b_;
  public:
    G96XFunctional();
    G96XFunctional(const RefKeyVal &);
    G96XFunctional(StateIn &);
    ~G96XFunctional();
    void save_data_state(StateOut &);

    int need_density_gradient();

    void point(const PointInputData&, PointOutputData&);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
