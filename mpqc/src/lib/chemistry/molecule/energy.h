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

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <math/optimize/function.h>
#include <math/optimize/conv.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/hess.h>


/** The MolecularEnergy abstract class inherits from the Function class.
It computes the energy of the molecule as a function of the geometry.  The
coordinate system used can be either internal or cartesian.  */
class MolecularEnergy: public Function {
#   define CLASSNAME MolecularEnergy
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension moldim_; // the number of cartesian variables
    RefMolecularCoor mc_;
    RefMolecule mol_;
    RefMolecularHessian hess_;
    RefMolecularHessian guesshess_;

    RefSCVector cartesian_gradient_;
    RefSymmSCMatrix cartesian_hessian_;
  protected:
    RefPointGroup initial_pg_;

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
  public:
    MolecularEnergy(const MolecularEnergy&);
    /** @memo The KeyVal constructor.
        \begin{description}

        \item[molecule] A Molecule object.  There is no default.

        \item[coor] A MolecularCoor object that describes the coordinates.
        If this is not given cartesian coordinates will be used.  For
        convenience, two keywords needed by the MolecularCoor object are
        automatically provided: natom3 and matrixkit.

        \item[value_accuracy] Sets the accuracy to which values are
        computed.  The default is 1.0e-6 atomic units.

        \item[gradient_accuracy] Sets the accuracy to which gradients are
        computed.  The default is 1.0e-6 atomic units.

        \item[hessian_accuracy] Sets the accuracy to which hessians are
        computed.  The default is 1.0e-4 atomic units.

        \item[print_molecule_when_changed] If true, then whenever the
        molecule's coordinates are updated they will be printed.  The
        default is true.

        \end{description}
    */
    MolecularEnergy(const RefKeyVal&);
    MolecularEnergy(StateIn&);
    ~MolecularEnergy();

    void save_data_state(StateOut&);

    MolecularEnergy & operator=(const MolecularEnergy&);
    
    /// A wrapper around value();
    virtual double energy();

    virtual RefMolecule molecule() const;
    virtual RefSCDimension moldim() const;
    
    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    /** If a molecule hessian object is given, it will be used to provide a
        hessian. */
    RefSymmSCMatrix hessian();
    int hessian_implemented() const;

    void set_x(const RefSCVector&);

    /// Return the cartesian coordinates.
    RefSCVector get_cartesian_x();
    /// Return the cartesian gradient.
    RefSCVector get_cartesian_gradient();
    /// Return the cartesian hessian.
    RefSymmSCMatrix get_cartesian_hessian();

    RefMolecularCoor molecularcoor() { return mc_; }

    /** Call this if you have changed the molecular symmetry of the
        molecule contained by this MolecularEnergy. */
    virtual void symmetry_changed();

    RefNonlinearTransform change_coordinates();
    
    /// Nicely print n x 3 data that are stored in a vector.
    void print_natom_3(const RefSCVector &,
                       const char *t=0, ostream&o=ExEnv::out()) const;
    void print_natom_3(double **, const char *t=0, ostream&o=ExEnv::out()) const;
    void print_natom_3(double *, const char *t=0, ostream&o=ExEnv::out()) const;

    virtual void print(ostream& = ExEnv::out()) const;
};
SavableState_REF_dec(MolecularEnergy);

class SumMolecularEnergy: public MolecularEnergy {
#   define CLASSNAME SumMolecularEnergy
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    RefMolecularEnergy *mole_;
    double *coef_;
    void compute();
  public:
    SumMolecularEnergy(const RefKeyVal &);
    SumMolecularEnergy(StateIn&);
    ~SumMolecularEnergy();

    void save_data_state(StateOut&);

    int value_implemented() const;
    int gradient_implemented() const;
    int hessian_implemented() const;

    void set_x(const RefSCVector&);
};
SavableState_REF_dec(SumMolecularEnergy);

/* The MolEnergyConvergence class derives from the Convergence class.  The
MolEnergyConvergence class allows the user to request that cartesian
coordinates be used in evaluating the convergence criteria.  This is
useful, since the internal coordinates can be somewhat arbitary.  If the
optimization is constrained, then the fixed internal coordinates will be
projected out of the cartesian gradients.  The input is similar to that for
Convergence class with the exception that giving none of the convergence
criteria keywords is the same as providing the following input to the
KeyVal constructor:

\begin{verbatim}
  conv<MolEnergyConvergence>: (
    max_disp = 1.0e-4
    max_grad = 1.0e-4
    graddisp = 1.0e-4
  )
\end{verbatim}

For MolEnergyConverence to work, the Function object given to the Optimizer
object must derive from MolecularEnergy.
*/
class MolEnergyConvergence: public Convergence {
#   define CLASSNAME MolEnergyConvergence
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefMolecularEnergy mole_;
    int cartesian_;

    void set_defaults();
  public:
    // Standard constructors and destructor.
    MolEnergyConvergence();
    MolEnergyConvergence(StateIn&);
    /** @memo The KeyVal constructor.

        In addition to the keywords read by Convergence, the following
        keywords are examined:

        \begin{description}

        \item[energy] The MolecularEnergy object.  This is required.

        \item[cartesian] If true, cartesian displacements and gradients
        will be compared to the convergence criteria.  The default is true.

        \end{description}

     */
    MolEnergyConvergence(const RefKeyVal&);
    virtual ~MolEnergyConvergence();

    void save_data_state(StateOut&);

    // Set the current gradient and position information.  These
    //will possibly grab the cartesian infomation if we have a
    //MolecularEnergy.
    void get_grad(const RefFunction &);
    void get_x(const RefFunction &);
    void set_nextx(const RefSCVector &);

    // Return nonzero if the optimization has converged.
    int converged();
};
SavableState_REF_dec(MolEnergyConvergence);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
