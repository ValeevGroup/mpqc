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
#include <util/group/thread.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/basis/extent.h>
#include <chemistry/qc/wfn/density.h>

namespace sc {

/** An abstract base class for integrating the electron density. */
class DenIntegrator: virtual public SavableState {
  protected:
    Ref<Wavefunction> wfn_;
//clj    Ref<ShellExtent> extent_;
    Ref<BatchElectronDensity> den_;

    Ref<ThreadGrp> threadgrp_;
    Ref<MessageGrp> messagegrp_;

    double value_;
    double accuracy_;

    double *alpha_vmat_;
    double *beta_vmat_;

//clj    double *alpha_dmat_;
//clj    double *beta_dmat_;
//clj    double *dmat_bound_;

    int spin_polarized_;

    int need_density_; // specialization must set to 1 if it needs density_
    double density_;
    int nbasis_;
    int nshell_;
    int n_integration_center_;
    int natom_;
    int compute_potential_integrals_; // 1 if potential integrals are needed

    int linear_scaling_;
    int use_dmat_bound_;

    void init_integration(const Ref<DenFunctional> &func,
                          const RefSymmSCMatrix& densa,
                          const RefSymmSCMatrix& densb,
                          double *nuclear_gradient);
    void done_integration();

    void init_object();
  public:
    /// Construct a new DenIntegrator.
    DenIntegrator();
    /// Construct a new DenIntegrator given the KeyVal input.
    DenIntegrator(const Ref<KeyVal> &);
    /// Construct a new DenIntegrator given the StateIn data.
    DenIntegrator(StateIn &);
    ~DenIntegrator();
    void save_data_state(StateOut &);

    /// Returns the wavefunction used for the integration.
    Ref<Wavefunction> wavefunction() const { return wfn_; }
    /// Returns the result of the integration.
    double value() const { return value_; }

    /// Sets the accuracy to use in the integration.
    void set_accuracy(double a);
    double get_accuracy(void) {return accuracy_; }
    /** Call with non zero if the potential integrals are to be computed.
        They can be returned with the vmat() member. */
    void set_compute_potential_integrals(int);
    /** Returns the alpha potential integrals. Stored as
        the lower triangular, row-major format. */
    const double *alpha_vmat() const { return alpha_vmat_; }
    /** Returns the beta potential integrals. Stored as
        the lower triangular, row-major format. */
    const double *beta_vmat() const { return beta_vmat_; }

    /** Called before integrate.  Does not need to be called again
        unless the geometry changes or done is called. */
    virtual void init(const Ref<Wavefunction> &);
    /// Must be called between calls to init.
    virtual void done();
    /** Performs the integration of the given functional using the given
        alpha and beta density matrices.  The nuclear derivative
        contribution is placed in nuclear_grad, if it is non-null. */
    virtual void integrate(const Ref<DenFunctional> &,
                           const RefSymmSCMatrix& densa =0,
                           const RefSymmSCMatrix& densb =0,
                           double *nuclear_grad = 0) = 0;
};


/** An abstract base class for computing grid weights. */
class IntegrationWeight: virtual public SavableState {

  protected:

    Ref<Molecule> mol_;
    double tol_;

    void fd_w(int icenter, SCVector3 &point, double *fd_grad_w);

  public:
    IntegrationWeight();
    IntegrationWeight(const Ref<KeyVal> &);
    IntegrationWeight(StateIn &);
    ~IntegrationWeight();
    void save_data_state(StateOut &);

    void test(int icenter, SCVector3 &point);
    void test();

    /// Initialize the integration weight object.
    virtual void init(const Ref<Molecule> &, double tolerance);
    /// Called when finished with the integration weight object.
    virtual void done();
    /** Computes the weight for a given center at a given point in space.
        Derivatives of the weigth with respect to nuclear coordinates are
        optionally returned in grad_w.  This must be called after init, but
        before done. It must also be thread-safe. */
    virtual double w(int center, SCVector3 &point, double *grad_w = 0) = 0;
};


/** Implements Becke's integration weight scheme. */
class BeckeIntegrationWeight: public IntegrationWeight {

    int n_integration_centers;
    SCVector3 *centers;
    double *atomic_radius;

    double **a_mat;
    double **oorab;

    void compute_grad_p(int gc, int ic, int wc, SCVector3&r, double p,
                           SCVector3&g);
    void compute_grad_nu(int gc, int bc, SCVector3 &point, SCVector3 &grad);

    double compute_t(int ic, int jc, SCVector3 &point);
    double compute_p(int icenter, SCVector3&point);

  public:
    BeckeIntegrationWeight();
    BeckeIntegrationWeight(const Ref<KeyVal> &);
    BeckeIntegrationWeight(StateIn &);
    ~BeckeIntegrationWeight();
    void save_data_state(StateOut &);

    void init(const Ref<Molecule> &, double tolerance);
    void done();
    double w(int center, SCVector3 &point, double *grad_w = 0);
};

/** An abstract base class for radial integrators. */
class RadialIntegrator: virtual public SavableState {
  protected:
    int nr_;
  public:
    RadialIntegrator();
    RadialIntegrator(const Ref<KeyVal> &);
    RadialIntegrator(StateIn &);
    ~RadialIntegrator();
    void save_data_state(StateOut &);

    /// The number of quadrature points (redundant since radial_value takes this as a parameter)
    virtual int nr() const = 0;
    /** returns the radius for the quadrature point ir.
	multiplier returns the quadrature weight which includes the Jacobian (r^2).
    */
    virtual double radial_value(int ir, int nr, double radii,
                                double &multiplier) = 0;
};


/** An abstract base class for angular integrators. */
class AngularIntegrator: virtual public SavableState{
  protected:
  public:
    AngularIntegrator();
    AngularIntegrator(const Ref<KeyVal> &);
    AngularIntegrator(StateIn &);
    ~AngularIntegrator();
    void save_data_state(StateOut &);

    virtual int nw(void) const = 0;
    virtual int num_angular_points(double r_value, int ir) = 0;
    virtual double angular_point_cartesian(int iangular, double r,
        SCVector3 &integration_point) const = 0;
};


/** An implementation of a radial integrator using the Euler-Maclaurin
    weights and grid points. */
class EulerMaclaurinRadialIntegrator: public RadialIntegrator {
  public:
    EulerMaclaurinRadialIntegrator();
    EulerMaclaurinRadialIntegrator(int i);
    /** Constructs a EulerMaclaurinRadialIntegrator from KeyVal input.
        The <tt>nr</tt> keyword gives the number of radial integration
        points.  The default is 75. */
    EulerMaclaurinRadialIntegrator(const Ref<KeyVal> &);
    EulerMaclaurinRadialIntegrator(StateIn &);
    ~EulerMaclaurinRadialIntegrator();
    void save_data_state(StateOut &);

    int nr() const;
    double radial_value(int ir, int nr, double radii, double &multiplier);

    void print(std::ostream & =ExEnv::out0()) const;
};

/** An implementation of a Lebedev angular integrator.  It uses code
    written by Dr. Dmitri N. Laikov.

    This can generate grids with the following numbers of points: 6, 14,
    26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 386, 434,
    482, 530, 590, 650, 698, 770, 830, 890, 974, 1046, 1118, 1202, 1274,
    1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222, 2354, 2450,
    2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470, 3590, 3722, 3890, 4010,
    4154, 4334, 4466, 4610, 4802, 4934, 5090, 5294, 5438, 5606, and 5810.

    V.I. Lebedev, and D.N. Laikov
    "A quadrature formula for the sphere of the 131st
     algebraic order of accuracy"
    Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   
    V.I. Lebedev
    "A quadrature formula for the sphere of 59th algebraic
     order of accuracy"
    Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
   
    V.I. Lebedev, and A.L. Skorokhodov
    "Quadrature formulas of orders 41, 47, and 53 for the sphere"
    Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
   
    V.I. Lebedev
    "Spherical quadrature formulas exact to orders 25-29"
    Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
   
    V.I. Lebedev
    "Quadratures on a sphere"
    Computational Mathematics and Mathematical Physics, Vol. 16,
    1976, pp. 10-24.
   
    V.I. Lebedev
    "Values of the nodes and weights of ninth to seventeenth
     order Gauss-Markov quadrature formulae invariant under the
     octahedron group with inversion"
    Computational Mathematics and Mathematical Physics, Vol. 15,
       1975, pp. 44-51.

 */
class LebedevLaikovIntegrator: public AngularIntegrator {
  protected:
    int npoint_;
    double *x_, *y_, *z_, *w_;
    
    void init(int n);
  public:
    LebedevLaikovIntegrator();
    /**
       Construct a LebedevLaikovIntegrator using the given KeyVal input.
       The <tt>n</tt> keyword gives the number of angular points.  The
       default is 302.
    */
    LebedevLaikovIntegrator(const Ref<KeyVal> &);
    LebedevLaikovIntegrator(StateIn &);
    LebedevLaikovIntegrator(int);
    ~LebedevLaikovIntegrator();
    void save_data_state(StateOut &);

    int nw(void) const;
    int num_angular_points(double r_value, int ir);
    double angular_point_cartesian(int iangular, double r,
                                   SCVector3 &integration_point) const;
    void print(std::ostream & =ExEnv::out0()) const;
};

/** An implementation of an angular integrator using the Gauss-Legendre
    weights and grid points. */
class GaussLegendreAngularIntegrator: public AngularIntegrator {
  protected:
    int ntheta_;
    int nphi_;
    int Ktheta_;
    int ntheta_r_;
    int nphi_r_;
    int Ktheta_r_;
    double *theta_quad_weights_;
    double *theta_quad_points_;

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
    int nw(void) const;
    double sin_theta(SCVector3 &point) const;
    void gauleg(double x1, double x2, int n);    
  public:
    GaussLegendreAngularIntegrator();
    /**
       Contract a GaussLegendreAngularIntegrator from KeyVal input.
       This class is for testing, the LebedevLaikovIntegrator
       is preferred for normal use.  The following parameters
       are read: <tt>ntheta</tt>, <tt>nphi</tt>, and <tt>Ktheta</tt>.
    */
    GaussLegendreAngularIntegrator(const Ref<KeyVal> &);
    GaussLegendreAngularIntegrator(StateIn &);
    ~GaussLegendreAngularIntegrator();
    void save_data_state(StateOut &);
    
    int num_angular_points(double r_value, int ir);
    double angular_point_cartesian(int iangular, double r,
        SCVector3 &integration_point) const;
    void print(std::ostream & =ExEnv::out0()) const;
};

/** An implementation of an integrator using any combination of
    a RadialIntegrator and an AngularIntegrator. */
class RadialAngularIntegrator: public DenIntegrator {
  private:
    int prune_grid_;
    double **Alpha_coeffs_;
    int gridtype_;
    int **nr_points_, *xcoarse_l_;
    int npruned_partitions_;
    double *grid_accuracy_;
    int dynamic_grids_;
    int natomic_rows_;
    int max_gridtype_;
  protected:
    Ref<IntegrationWeight> weight_;
    Ref<RadialIntegrator> radial_user_;
    Ref<AngularIntegrator> angular_user_;
    Ref<AngularIntegrator> ***angular_grid_;
    Ref<RadialIntegrator> **radial_grid_;
  public:
    RadialAngularIntegrator();
    /** Construct a RadialAngularIntegrator from KeyVal input.

        The accepted keyword are listed below.  The most important keyword
        is <tt>grid</tt>.  The <tt>dynamic</tt> and <tt>prune_grid</tt>
        options may be of occassional interest.
        <dl>

        <dt><tt>grid</tt><dd>Specifies the fineness of the grid.  Possible
        values are <tt>xcoarse</tt>, <tt>coarse</tt>, <tt>medium</tt>,
        <tt>fine</tt>, <tt>xfine</tt>, and <tt>ultrafine</tt>, in order of
        increasing accuracy and cost.  The default is <tt>fine</tt>.

        <dt><tt>dynamic</tt><dd>This gives a boolean value that, if true,
        will cause the grids to start out coarse, and approach the
        requested <tt>grid</tt> value as more accuracy is required, when
        the calculation is close to convergence.  The default is true.

        <dt><tt>prune_grid</tt><dd>This gives a boolean value that, if
        true, will cause more course angular grids to be used near
        nuclei.  The default is true.  When this is true, further control
        over pruning can be obtained with the <tt>angular_points</tt>
        and <tt>alpha_coeffs</tt> keywords.

        <dt><tt>radial</tt><dd>Specifies the RadialIntegrator object.  If
        this is given, then specifying the <tt>grid</tt> and
        <tt>dynamic</tt> keywords will not affect the radial grid.  The
        default is controlled by other options, but is always one of
        several EulerMaclaurinRadialIntegrator objects.

        <dt><tt>angular</tt><dd>Specifies the AngularIntegrator object.  If
        this is given, then specifying the <tt>grid</tt>,
        <tt>prune_grid</tt>, and <tt>dynamic</tt> keywords will not affect
        the angular grid.  The default is controlled by other options,
        but is always one of several LebedevLaikovIntegrator objects.

        <dt><tt>weight</tt><dd>Specifies the IntegrationWeight object.
        The default is BeckeIntegrationWeight.

        </dl>
     */
    RadialAngularIntegrator(const Ref<KeyVal> &);
    RadialAngularIntegrator(StateIn &);
    ~RadialAngularIntegrator();
    void save_data_state(StateOut &);

    void integrate(const Ref<DenFunctional> &,
                   const RefSymmSCMatrix& densa =0,
                   const RefSymmSCMatrix& densb =0,
                   double *nuclear_gradient = 0);
    void print(std::ostream & =ExEnv::out0()) const;
    AngularIntegrator *get_angular_grid(double radius, double atomic_radius,
                                        int charge, int deriv_order);
    RadialIntegrator *get_radial_grid(int charge, int deriv_order);
    void init_default_grids(void);
    int angular_grid_offset(int i);
    void set_grids(void);
    int get_atomic_row(int i);
    void init_parameters(void);
    void init_parameters(const Ref<KeyVal>& keyval);
    void init_pruning_coefficients(const Ref<KeyVal>& keyval);
    void init_pruning_coefficients(void);
    void init_alpha_coefficients(void);
    int select_dynamic_grid(void);
    Ref<IntegrationWeight> weight() { return weight_; }
};

}
    
#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
