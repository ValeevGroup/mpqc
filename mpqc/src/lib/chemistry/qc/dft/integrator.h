//
// integrator.h --- definition of the electron density integrator
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

#ifndef _chemistry_qc_dft_integrator_h
#define _chemistry_qc_dft_integrator_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/basis/extent.h>

class DenIntegrator: virtual public SavableState {
#   define CLASSNAME DenIntegrator
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefWavefunction wfn_;
    RefShellExtent extent_;

    double value_;
    double accuracy_;

    int spin_polarized_;

    int ncontrib_;
    int *contrib_;
    int ncontrib_bf_;
    int *contrib_bf_;
    double *bs_values_;
    double *bsg_values_;
    double *bsh_values_;
    double *alpha_dmat_;
    double *beta_dmat_;
    double *dmat_bound_;
    double *alpha_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    double *beta_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    int need_density_; // specialization must set to 1 if it needs density_
    double density_;
    int nbasis_;
    int nshell_;
    int natom_;
    int compute_potential_integrals_; // 1 if potential integrals are needed

    int need_gradient_;
    int need_hessian_;

    int linear_scaling_;
    int use_dmat_bound_;

    void get_density(double *dmat, PointInputData::SpinData &d);
    void init_integration(const RefDenFunctional &func,
                          const RefSymmSCMatrix& densa,
                          const RefSymmSCMatrix& densb,
                          double *nuclear_gradient);
    void done_integration();
    double do_point(int acenter, const SCVector3 &r, const RefDenFunctional &,
                    double weight, double multiplier, double *nuclear_gradient,
                    double *f_gradient, double *w_gradient);
  public:
    DenIntegrator();
    DenIntegrator(const RefKeyVal &);
    DenIntegrator(StateIn &);
    ~DenIntegrator();
    void save_data_state(StateOut &);

    RefWavefunction wavefunction() const { return wfn_; }
    double value() const { return value_; }

    void set_accuracy(double a) { accuracy_ = a; }

    // Call with non zero if the potential integrals are to be computed.
    // They can be returned with the vmat() member.
    void set_compute_potential_integrals(int);
    const double *alpha_vmat() const { return alpha_vmat_; }
    const double *beta_vmat() const { return beta_vmat_; }

    /** Called before integrate.  Does not need to be called again
        unless the geometry changes or done is called. */
    virtual void init(const RefWavefunction &);
    /// Must be called between calls to init.
    virtual void done();
    virtual void integrate(const RefDenFunctional &,
                           const RefSymmSCMatrix& densa =0,
                           const RefSymmSCMatrix& densb =0,
                           double *nuclear_grad = 0) = 0;
};
SavableState_REF_dec(DenIntegrator);

class IntegrationWeight: virtual public SavableState {
#   define CLASSNAME IntegrationWeight
#   include <util/state/stated.h>
#   include <util/class/classda.h>

  protected:

    RefMolecule mol_;
    double tol_;

    void fd_w(int icenter, SCVector3 &point, double *fd_grad_w);

  public:
    IntegrationWeight();
    IntegrationWeight(const RefKeyVal &);
    IntegrationWeight(StateIn &);
    ~IntegrationWeight();
    void save_data_state(StateOut &);

    void test(int icenter, SCVector3 &point);
    void test();

    virtual void init(const RefMolecule &, double tolerance);
    virtual void done();
    virtual double w(int center, SCVector3 &point, double *grad_w = 0) = 0;
};
SavableState_REF_dec(IntegrationWeight);

class BeckeIntegrationWeight: public IntegrationWeight {
#   define CLASSNAME BeckeIntegrationWeight
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>

    int ncenters;
    SCVector3 *centers;
    double *bragg_radius;

    double **a_mat;
    double **oorab;

    void compute_grad_p(int gc, int ic, int wc, SCVector3&r, double p,
                           SCVector3&g);
    void compute_grad_nu(int gc, int bc, SCVector3 &point, SCVector3 &grad);

    double compute_t(int ic, int jc, SCVector3 &point);
    double compute_p(int icenter, SCVector3&point);

  public:
    BeckeIntegrationWeight();
    BeckeIntegrationWeight(const RefKeyVal &);
    BeckeIntegrationWeight(StateIn &);
    ~BeckeIntegrationWeight();
    void save_data_state(StateOut &);

    void init(const RefMolecule &, double tolerance);
    void done();
    double w(int center, SCVector3 &point, double *grad_w = 0);
};

class RadialIntegrator: virtual public SavableState{
#   define CLASSNAME RadialIntegrator
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int nr_;
  public:
    RadialIntegrator();
    RadialIntegrator(const RefKeyVal &);
    RadialIntegrator(StateIn &);
    ~RadialIntegrator();
    void save_data_state(StateOut &);

    void set_nr(int i);
    int get_nr(void) const;
    virtual double radial_value(int ir, int nr, double radii) = 0;
    virtual double radial_multiplier(int nr) = 0;
    virtual double get_dr_dq(void) const = 0;
    virtual double get_dr_dqr2(void) const = 0;
    virtual void set_dr_dq(double i) = 0;
    virtual void set_dr_dqr2(double i) = 0;
    void print(ostream & =cout) const;
};
SavableState_REF_dec(RadialIntegrator);

class AngularIntegrator: virtual public SavableState{
#   define CLASSNAME AngularIntegrator
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
  public:
    AngularIntegrator();
    AngularIntegrator(const RefKeyVal &);
    AngularIntegrator(StateIn &);
    ~AngularIntegrator();
    void save_data_state(StateOut &);

    virtual int num_angular_points(double r_value, int ir) = 0;
    virtual void angular_weights(void) = 0;
    virtual double angular_point_cartesian(int iangular, SCVector3 &point,
        SCVector3 &integration_point) const = 0;
    virtual void print(ostream & =cout) const = 0;
};
SavableState_REF_dec(AngularIntegrator);

class EulerMaclaurinRadialIntegrator: public RadialIntegrator {
#   define CLASSNAME EulerMaclaurinRadialIntegrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double dr_dq_;
    double dr_dqr2_;
  public:
    EulerMaclaurinRadialIntegrator();
    EulerMaclaurinRadialIntegrator(const RefKeyVal &);
    EulerMaclaurinRadialIntegrator(StateIn &);
    ~EulerMaclaurinRadialIntegrator();
    void save_data_state(StateOut &);

    double radial_value(int ir, int nr, double radii);
    double radial_multiplier(int nr);
    double get_dr_dq(void) const;
    void set_dr_dq(double i);
    double get_dr_dqr2(void) const;
    void set_dr_dqr2(double i);
};

class LebedevAngularIntegrator: public AngularIntegrator {
#   define CLASSNAME LebedevAngularIntegrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int norder_;
    int npoints_;
    int N1_, N2_, N3_;
    double *x_, *y_, *z_;
    double *lebedev_weights_;
    int point_count_;
  public:
    LebedevAngularIntegrator();
    LebedevAngularIntegrator(const RefKeyVal &);
    LebedevAngularIntegrator(StateIn &);
    ~LebedevAngularIntegrator();
    void save_data_state(StateOut &);

    int get_norder(void) const;
    void set_norder(int i);
    int get_npoints(void) const;
    void set_npoints(int i);
    int get_N1(void) const;
    void set_N1(int i);
    int get_N2(void) const;
    void set_N2(int i);
    int get_N3(void) const;
    void set_N3(int i);
    int get_point_count(void) const;
    void set_point_count(int i);
    double angular_point_cartesian(int iangular, SCVector3 &point,
                                   SCVector3 &integration_point) const;
    int num_angular_points(double r_value, int ir);
    void angular_weights(void);
    void build_grid(void);
    void generate_points(double weights[], int N, int nsets, double  u[], double v[], double w[]);
    void expand(double array[], int offset, double weight);
    void print(ostream & =cout) const;
};

class GaussLegendreAngularIntegrator: public AngularIntegrator {
#   define CLASSNAME GaussLegendreAngularIntegrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int ntheta_;
    int nphi_;
    int Ktheta_;
    int ntheta_r_;
    int nphi_r_;
    int Ktheta_r_;
    double *theta_quad_weights_;
    double *theta_quad_points_;
  public:
    GaussLegendreAngularIntegrator();
    GaussLegendreAngularIntegrator(const RefKeyVal &);
    GaussLegendreAngularIntegrator(StateIn &);
    ~GaussLegendreAngularIntegrator();
    void save_data_state(StateOut &);
    
    int get_ntheta(void) const;
    void set_ntheta(int i);
    int get_nphi(void) const;
    void set_nphi(int i);
    int get_Ktheta(void) const;
    void set_Ktheta(int i);
    int get_ntheta_r(void) const;
    void set_ntheta_r(int i);
    int get_nphi_r(void) const;
    void set_nphi_r(int i);
    int get_Ktheta_r(void) const;
    void set_Ktheta_r(int i);
    int num_angular_points(double r_value, int ir);
    void angular_weights(void);
    double angular_point_cartesian(int iangular, SCVector3 &point,
        SCVector3 &integration_point) const;
    double sin_theta(SCVector3 &point) const;
    void gauleg(double x1, double x2, int n);    
    void print(ostream & =cout) const;
};

class RadialAngularIntegrator: public DenIntegrator {
#   define CLASSNAME RadialAngularIntegrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefRadialIntegrator RadInt_;
    RefAngularIntegrator AngInt_;
    RefIntegrationWeight weight_;
  public:
    RadialAngularIntegrator();
    RadialAngularIntegrator(const RefKeyVal &);
    RadialAngularIntegrator(StateIn &);
    ~RadialAngularIntegrator();
    void save_data_state(StateOut &);

    void integrate(const RefDenFunctional &,
                   const RefSymmSCMatrix& densa =0,
                   const RefSymmSCMatrix& densb =0,
                   double *nuclear_gradient = 0);

    void print(ostream & =cout) const;
};
    
// Based on C.W. Murray, et al. Mol. Phys. 78, No. 4, 997-1014, 1993
class Murray93Integrator: public DenIntegrator {
#   define CLASSNAME Murray93Integrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int nr_;
    int ntheta_;
    int nphi_;
    int Ktheta_;

    RefIntegrationWeight weight_;
    
  public:
    Murray93Integrator();
    Murray93Integrator(const RefKeyVal &);
    Murray93Integrator(StateIn &);
    ~Murray93Integrator();
    void save_data_state(StateOut &);

    void integrate(const RefDenFunctional &,
                   const RefSymmSCMatrix& densa =0,
                   const RefSymmSCMatrix& densb =0,
                   double *nuclear_gradient = 0);

    void print(ostream & =cout) const;
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
