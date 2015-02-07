//
// energy.h
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

#ifndef _chemistry_molecule_energy_h
#define _chemistry_molecule_energy_h

#include <iostream>

#include <math/optimize/function.h>
#include <math/optimize/conv.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/deriv.h>
#ifdef MPQC_NEW_FEATURES
#include <util/misc/xml.h>
#endif

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

/** The MolecularEnergy abstract class inherits from the Function class.
It computes the energy of the molecule as a function of the geometry.  The
coordinate system used can be either internal or cartesian.  */
class MolecularEnergy: public Function
#ifdef MPQC_NEW_FEATURES
, virtual public DescribedXMLWritable
#endif
{
  private:
    RefSCDimension moldim_; // the number of cartesian variables
    Ref<MolecularCoor> mc_;
    Ref<Molecule> mol_;
    /** it seems to be a bad idea to have this here -- in order to initialize hess I may need to
        call virtual functions of this, which are not available in the constructor
     */
    Ref<MolecularHessian> hess_;
    Ref<MolecularHessian> guesshess_;
    Ref<MolecularGradient> grad_;

    RefSCVector cartesian_gradient_;
    RefSymmSCMatrix cartesian_hessian_;

    RefSCVector efield_;  //< electric field vector

    /// Whether to do intermediate checkpointing of this object
    bool ckpt_;
    /// Name of the file into which to checkpoint this object
    std::string ckpt_file_;
    /// How often this object should be checkpointed (only matters in iterative methods)
    int ckpt_freq_;

  protected:
    Ref<PointGroup> initial_pg_;

    void failure(const char *);

    /// This is just a wrapper around set_value().
    virtual void set_energy(double);

    /** These are passed gradients and hessian in cartesian coordinates.
        The gradient and hessian in internal coordinates are computed. */
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    void x_to_molecule();
    void molecule_to_x();

    int print_molecule_when_changed_;

    /** overload this in classes that support computations in nonzero electric field
        the default is to not support external electric fields. */
    virtual bool nonzero_efield_supported() const;

    /** must overload this in a derived class if analytic gradient can be computed
     * @return true (analytic gradient is available) or false (analytic gradient is not available, default)
     */
    virtual bool analytic_gradient_implemented() const;
    /** must overload this in a derived class if analytic hessian can be computed
     * @return true (analytic hessian is available) or false (analytic hessian is not available, default)
     */
    virtual bool analytic_hessian_implemented() const;


  public:
    MolecularEnergy(const MolecularEnergy&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>molecule</tt><dd> A Molecule object.  There is no default.

        <dt><tt>coor</tt><dd> A MolecularCoor object that describes the
        coordinates.  If this is not given cartesian coordinates will be
        used.  For convenience, two keywords needed by the MolecularCoor
        object are automatically provided: <tt>natom3</tt> and <tt>matrixkit</tt>.

        <dt><tt>value_accuracy</tt><dd> Sets the accuracy to which values
        are computed.  The default is 1.0e-6 atomic units.

        <dt><tt>gradient_accuracy</tt><dd> Sets the accuracy to which
        gradients are computed.  The default is 1.0e-6 atomic units.

        <dt><tt>hessian_accuracy</tt><dd> Sets the accuracy to which
        hessians are computed.  The default is 1.0e-4 atomic units.

        <dt><tt>hessian</tt><dd>Specifies a MolecularHessian object that is
        used to compute the hessian. This keyword may only need to be specified.
        if "exact" hessian is needed but this MolecularEnergy
        specialization does not provide a hessian of its own.

        <dt><tt>guess_hessian</tt><dd>Specifies a MolecularHessian object
        that is used to compute a guess hessian.  Guess hessians are used
        to improve the rate of convergence of optimizations.  If this
        keyword is not specified, and a MolecularCoor object is given by
        <tt>coor</tt>, then the guess hessian is obtained from the
        MolecularCoor object.  If neither this nor <tt>coor</tt> are given,
        then Function::guess_hessian is used, which returns a unit matrix.

        <dt><tt>electric_field</tt><dd> This 3-element vector specifies the Cartesian
        components of an external uniform electric field, in a.u. The default value, <tt>[0 0 0]</tt>,
        causes computations in absence of an electric field. Not all MolecularEnergy objects
        will support computations in presence of an electric field -- use
        MolecularEnergy::nonzero_efield_supported() to query objects about this capability.

        <dt><tt>print_molecule_when_changed</tt><dd> If true, then whenever
        the molecule's coordinates are updated they will be printed.  The
        default is true.

	<dt><tt>checkpoint</tt><dd> If true, then this object will be
	checkpointed during its evaluation. Not all implementations
	of <tt>MolecularEnergy</tt> support checkpointing.
	The default is false.

	<dt><tt>checkpoint_file</tt><dd> Specifies the name of the file
	into which this object will be checkpointed. Default is
	"<inpubasename>.ckpt", where "<inputbasename>" is the name of the input
	file without ".in".

	<dt><tt>checkpoint_freq</tt><dd> Specifies how often this object to
	be checkpointed. Only matters for objects which are computed
	iteratively. Default is 1.
	</dl> */
    MolecularEnergy(const Ref<KeyVal>&);
    MolecularEnergy(StateIn&);
    ~MolecularEnergy();

    void save_data_state(StateOut&);


#ifdef MPQC_NEW_FEATURES
    virtual boost::property_tree::ptree& write_xml(boost::property_tree::ptree& parent, const XMLWriter& writer);
#endif

    /// Set up checkpointing
    void set_checkpoint();
    void set_checkpoint_file(const char*);
    void set_checkpoint_freq(int freq);
    /// Check if need to checkpoint
    bool if_to_checkpoint() const;
    const char* checkpoint_file() const;
    int checkpoint_freq() const;

    MolecularEnergy & operator=(const MolecularEnergy&);

    /// A wrapper around value();
    virtual double energy();

    virtual Ref<Molecule> molecule() const;
    virtual RefSCDimension moldim() const;

    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    /** Reports whether gradient is implemented either analytically or using MolecularGradient object.
     * I don't see a need to reimplement this in a derived class
     * @return 0 (gradient cannot be computed) or 1 (gradient can be computed)
     */
    int gradient_implemented() const;
    /** Reports whether hessian is implemented either analytically or using MolecularHessian object.
     * I don't see a need to reimplement this in a derived class
     * @return 0 (hessian cannot be computed) or 1 (hessian can be computed)
     */
    int hessian_implemented() const;

    /**
     * These functions overload their Function counterparts.
     * If hessian/gradient objects are provided, these functions will convey desired accuracy to them.
     */
    //@{
    void set_desired_gradient_accuracy(double acc);
    void set_desired_hessian_accuracy(double acc);
    //@}

    /// Use this function to provide MolecularHessian object
    /// that will be used to compute hessian. You can call this function with null pointer to restore the state
    /// to the original state.
    void set_molhess(const Ref<MolecularHessian>& molhess);
    const Ref<MolecularHessian>& molhess() const;
    /// Will throw if hessian_implemented() returns 0
    RefSymmSCMatrix hessian();
    /// Use this function to provide MolecularGradient object
    /// that will be used to compute gradient. You can call this function with null pointer to restore the state
    /// to the original state.
    void set_molgrad(const Ref<MolecularGradient>& molgrad);
    const Ref<MolecularGradient>& molgrad() const;
    /// Will throw if gradient_implemented() returns 0
    RefSCVector gradient();

    void set_x(const RefSCVector&);

    /// Return the cartesian coordinates.
    RefSCVector get_cartesian_x();
    /// Return the cartesian gradient.
    RefSCVector get_cartesian_gradient();
    /// Return the cartesian hessian.
    RefSymmSCMatrix get_cartesian_hessian();

    Ref<MolecularCoor> molecularcoor() { return mc_; }

    /** Call this if you have changed the molecular symmetry of the
        molecule contained by this MolecularEnergy. */
    virtual void symmetry_changed();

    Ref<NonlinearTransform> change_coordinates();

    /** This function purges any caches of data in MolecularEnergy. It is useful with MolecularEnergy objects
        that keep state when obsolete() is called (for example, it makes sense for SCF to keep its old eigenvector
        and reuse it as a guess when geometry changes). The default implementation does nothing
        and must be overloaded in classes which need it */
    virtual void purge();

    /// returns the electric field vector
    const RefSCVector& electric_field() const { return efield_; }

    /// Nicely print n x 3 data that are stored in a vector.
    void print_natom_3(const RefSCVector &,
                       const char *t=0, std::ostream&o=ExEnv::out0()) const;
    void print_natom_3(double **, const char *t=0, std::ostream&o=ExEnv::out0()) const;
    void print_natom_3(double *, const char *t=0, std::ostream&o=ExEnv::out0()) const;

    virtual void print(std::ostream& = ExEnv::out0()) const;
};

/// linear combination of MolecularEnergy objects
class SumMolecularEnergy: public MolecularEnergy {
  protected:
    int value_implemented() const;
    bool analytic_gradient_implemented() const;
    bool analytic_hessian_implemented() const;
    int n_;
    Ref<MolecularEnergy> *mole_;
    double *coef_;
    void compute();
  public:
    SumMolecularEnergy(const Ref<KeyVal> &);
    SumMolecularEnergy(StateIn&);
    ~SumMolecularEnergy();

    void save_data_state(StateOut&);

    void set_x(const RefSCVector&);

    void purge();
};


/* The MolEnergyConvergence class derives from the Convergence class.  The
MolEnergyConvergence class allows the user to request that cartesian
coordinates be used in evaluating the convergence criteria.  This is
useful, since the internal coordinates can be somewhat arbitary.  If the
optimization is constrained, then the fixed internal coordinates will be
projected out of the cartesian gradients.  The input is similar to that for
Convergence class with the exception that giving none of the convergence
criteria keywords is the same as providing the following input to the
KeyVal constructor:

<pre>
  conv<MolEnergyConvergence>: (
    max_disp = 1.0e-4
    max_grad = 1.0e-4
    graddisp = 1.0e-4
  )
</pre>

For MolEnergyConverence to work, the Function object given to the Optimizer
object must derive from MolecularEnergy.
*/
class MolEnergyConvergence: public Convergence {
  protected:
    Ref<MolecularEnergy> mole_;
    int cartesian_;

    void set_defaults();
  public:
    // Standard constructors and destructor.
    MolEnergyConvergence();
    MolEnergyConvergence(StateIn&);
    /** The KeyVal constructor.

        In addition to the keywords read by Convergence, the following
        keywords are examined:

        <dl>

        <dt><tt>energy</tt><dd> The MolecularEnergy object.  This is
        required.

        <dt><tt>cartesian</tt><dd> If true, cartesian displacements and
        gradients will be compared to the convergence criteria.  The
        default is true.

        </dl>

     */
    MolEnergyConvergence(const Ref<KeyVal>&);
    virtual ~MolEnergyConvergence();

    void save_data_state(StateOut&);

    // Set the current gradient and position information.  These
    //will possibly grab the cartesian infomation if we have a
    //MolecularEnergy.
    void get_grad(const Ref<Function> &);
    void get_x(const Ref<Function> &);
    void set_nextx(const RefSCVector &);

    // Return nonzero if the optimization has converged.
    int converged();

    void print(std::ostream& = ExEnv::out0()) const;
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
