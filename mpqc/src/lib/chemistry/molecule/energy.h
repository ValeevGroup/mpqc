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

//. The \clsnm{MolecularEnergy} abstract class is a \clsnmref{Function}
//specialization that can optionally transform the cartesian molecule
//geometry, gradient, and hessian into an internal coordinate system.
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

    //. This is just a wrapper around set_value().
    virtual void set_energy(double);

    //. These are passed gradients and hessian in cartesian coordinates.
    //The gradient and hessian in internal coordinates are computed.
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    void x_to_molecule();
    void molecule_to_x();

    int print_molecule_when_changed_;
  public:
    MolecularEnergy(const MolecularEnergy&);
    MolecularEnergy(const RefKeyVal&);
    MolecularEnergy(StateIn&);
    ~MolecularEnergy();

    void save_data_state(StateOut&);

    MolecularEnergy & operator=(const MolecularEnergy&);
    
    //. A wrapper around value();
    virtual double energy();

    virtual RefMolecule molecule();
    virtual RefSCDimension moldim();
    
    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    //. If a molecule hessian object is given, it will be used
    //to provide a hessian.
    RefSymmSCMatrix hessian();
    int hessian_implemented();

    void set_x(const RefSCVector&);

    //. Return the cartesian coordinates.
    RefSCVector get_cartesian_x();
    //. Return the cartesian gradient.
    RefSCVector get_cartesian_gradient();
    //. Return the cartesian hessian.
    RefSymmSCMatrix get_cartesian_hessian();

    RefMolecularCoor molecularcoor() { return mc_; }

    //. Call this if you have changed the molecular symmetry of the
    // molecule contained by this MolecularEnergy
    virtual void symmetry_changed();

    RefNonlinearTransform change_coordinates();
    
    //. Nicely print n x 3 data that are stored in a vector.
    void print_natom_3(const RefSCVector &, const char *t=0, ostream&o=cout);

    virtual void print(ostream& = cout);
};
SavableState_REF_dec(MolecularEnergy);

class SumMolecularEnergy: public MolecularEnergy {
#   define CLASSNAME SumMolecularEnergy
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
    RefMolecularEnergy *mole_;
    double *coef_;
  protected:
    void compute();
  public:
    SumMolecularEnergy(const RefKeyVal &);
    SumMolecularEnergy(StateIn&);
    ~SumMolecularEnergy();

    void save_data_state(StateOut&);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();

    void set_x(const RefSCVector&);
};
SavableState_REF_dec(SumMolecularEnergy);

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
