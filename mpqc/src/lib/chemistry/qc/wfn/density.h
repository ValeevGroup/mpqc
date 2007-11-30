//
// density.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_wfn_density_h
#define _chemistry_qc_wfn_density_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/volume.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/basis/extent.h>
#include <chemistry/molecule/molrender.h>
#include <chemistry/molecule/molecule.h>
#include <math/mmisc/grid.h>

namespace sc {

/** This is a Volume that computes the electron density.  It
    can be used to generate isodensity surfaces. */
class ElectronDensity: public Volume {
  protected:
    Ref<Wavefunction> wfn_;
    virtual void compute();
  public:
    ElectronDensity(const Ref<KeyVal>&);
    ElectronDensity(const Ref<Wavefunction>&);
    ~ElectronDensity();
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2);
};

/** This a more highly optimized than ElectronDensity since
    everything is precomputed.  However, it cannot be used
    if the density and/or geometry might change between
    computations of the density or bounding box, unless the
    obsolete member is called. */
class BatchElectronDensity: public Volume {
    void zero_pointers();
  protected:
    // wfn_ might be null in which case basis_ and integral_
    // must be initialized and set_densities must be used.
    Ref<Wavefunction> wfn_;

    Ref<GaussianBasisSet> basis_;
    Ref<Integral> integral_;

    bool initialized_;

    // shared between threads
    double *alpha_dmat_;
    double *beta_dmat_;
    double *dmat_bound_;
    ShellExtent *extent_;

    // private data
    GaussianBasisSet::ValueData *valdat_;
    int ncontrib_;
    int *contrib_;
    int ncontrib_bf_;
    int *contrib_bf_;
    double *bs_values_;
    double *bsg_values_;
    double *bsh_values_;

    int nshell_;
    int nbasis_;
    bool spin_polarized_;
    int linear_scaling_;
    int use_dmat_bound_;

    bool need_hessian_, need_gradient_;
    bool need_basis_hessian_, need_basis_gradient_;

    bool using_shared_data_;

    double accuracy_;
    virtual void init_common_data();
    // this must be called after common data is initialized,
    // either with init_common_data or by copying
    virtual void init_scratch_data();
    void compute_basis_values(const SCVector3&r);
    void compute_spin_density(const double *restrictxx dmat,
                              double *restrictxx rho,
                              double *restrictxx grad,
                              double *restrictxx hess);

    virtual void compute();
  public:

    /** This gives the elements of gradient arrays. */
    enum {X=0, Y=1, Z=2};
    /** This gives the elements of hessian arrays. */
    enum {XX=0, YX=1, YY=2, ZX=3, ZY=4, ZZ=5};

    BatchElectronDensity(const Ref<KeyVal>&);
    BatchElectronDensity(const Ref<GaussianBasisSet> &basis,
                         const Ref<Integral> &integral,
                         double accuracy=DBL_EPSILON);
    BatchElectronDensity(const Ref<Wavefunction> &wfn,
                         double accuracy=DBL_EPSILON);
    /** This will construct copies of this.  If reference_parent_data is
        true, then data that do not change, such as the density matrices
        and shell extent, are referenced rather than copied.  In this case,
        the original object that allocated this items must be valid while
        copied objects are used to compute densities.  Also d must have
        already been intialized and the resulting copy is already initialized
        (and cannot be reinitialized).

        If reference_parent_data is false, then init must be called on this
        object before it is used.  */
    BatchElectronDensity(const Ref<BatchElectronDensity>& d,
                         bool reference_parent_data=false);
    ~BatchElectronDensity();
    /** Returns the bounding box. */
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2);

    /** This will cause all stratch storage to be released. */
    void clear();

    /** This is a alternate to the Volume interface that avoids some of the
        overhead of that interface. */
    void compute_density(const SCVector3 &r,
                         double *alpha_density,
                         double *alpha_density_grad,
                         double *alpha_density_hessian,
                         double *beta_density,
                         double *beta_density_grad,
                         double *beta_density_hessian);

    /** This is called to finish initialization of the object.  It must not
        be called with objects created in a way that they share parent
        data; those objects are initialized when they are constructed. This
        member is usually called automatically; however, for objects
        which will be used to provide other BatchElectronDensity objects'
        data, this must member be called before the other
        BatchElectronDensity objects are constructed.  If
        initialize_density_matrices is false, then the density matrices
        will be allocated, but not filled in.  They must be later filled in
        with set_densities. */
    virtual void init();

    /** This will fill in the internal copies of the density matrices with
        new values.  aden is the alpha density matrix and bden is the beta
        density matrix.  If bden is not null and does not reference the
        same matrix that is referenced by bden, then a spin polarized
        computation of the density is performed. */
    virtual void set_densities(const RefSymmSCMatrix &aden,
                               const RefSymmSCMatrix &bden = 0);

    /** Use the given Wavefunction object to set the electron density
        matrices. */
    virtual void set_densities(const Ref<Wavefunction> &wfn);

    /** Turn linear scaling algorithm on/off. The effect of this will be
        delayed until the next time init() is called. */
    void set_linear_scaling(bool b) { linear_scaling_ = b; }

    /** Sets the accuracy.  */
    void set_accuracy(double a) { accuracy_ = a; }

    /** Turn use of density matrix bounds on/off. */
    void set_use_dmat_bound(bool b) { use_dmat_bound_ = b; }

    /** Return true if the densities are spin polarized.  Set_densities
        must have been called, at least implicitly. */
    bool spin_polarized() const { return spin_polarized_; }

    /** @name DFT Support Members.
        These return some of the internal data, some of which is
        only value after a density has been computed.  This data
        is needed by the density functional theory code. */
    //@{
    /** Return the alpha density matrix. */
    double *alpha_density_matrix() { return alpha_dmat_; }
    /** Return the beta density matrix. */
    double *beta_density_matrix()
        { return (spin_polarized_?beta_dmat_:alpha_dmat_); }
    int ncontrib() { return ncontrib_; }
    int *contrib() { return contrib_; }
    int ncontrib_bf() { return ncontrib_bf_; }
    int *contrib_bf() { return contrib_bf_; }
    double *bs_values() { return bs_values_; }
    double *bsg_values() { return bsg_values_; }
    double *bsh_values() { return bsh_values_; }
    /** To ensure that that the basis functions gradients are computed,
        use this. */
    void set_need_basis_gradient(bool b) { need_basis_gradient_ = b; }
    void set_need_basis_hessian(bool b) { need_basis_hessian_ = b; }
    //@}
};

/** The WriteElectronDensity class writes the electron density at user defined
    grid points to the standard output or to a separate file. */
class WriteElectronDensity: public WriteGrid {
  private:
    double df_alpha(double alpha, double beta);
    double df_beta(double alpha, double beta);
    double df_sum(double alpha, double beta);
    double df_spin(double alpha, double beta);
  protected:
    Ref<Wavefunction> wfn_;
    Ref<BatchElectronDensity> bed_;
    double accuracy_;
    char *type_;
    double (WriteElectronDensity::*density_function_)(double, double);
    
    void initialize();
    void label(char* buffer);
    Ref<Molecule> get_molecule();
    double calculate_value(SCVector3 point);
  public:
    /** The KeyVal constructor
        
        <dl>

        <dt><tt>wfn</tt></dt><dd> The Wavefunction of which the density is
        calculated. There is no default for this option.</dd>

        <dt><tt>type</tt></dt><dd> Four types of densities can be written to
        the output: 'sum', 'alpha', 'beta' or 'spin'. The default is 'sum'.</dd>

        <dt><tt>accuracy</tt></dt><dd>The accuracy to which the density is
        calculated. The default is maximum accuracy.</dd>

        </dl> */
    WriteElectronDensity(const Ref<KeyVal> &);
};

class DensityColorizer: public MoleculeColorizer {
  protected:
    Ref<Wavefunction> wfn_;
    double scale_;
    double reference_;
    int have_scale_;
    int have_reference_;
  public:
    DensityColorizer(const Ref<KeyVal>&);
    ~DensityColorizer();

    void colorize(const Ref<RenderedPolygons> &);
};

class GradDensityColorizer: public MoleculeColorizer {
  protected:
    Ref<Wavefunction> wfn_;
    double scale_;
    double reference_;
    int have_scale_;
    int have_reference_;
  public:
    GradDensityColorizer(const Ref<KeyVal>&);
    ~GradDensityColorizer();

    void colorize(const Ref<RenderedPolygons> &);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
