//
// energy.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

class MolecularEnergy: public Function {
#   define CLASSNAME MolecularEnergy
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension moldim_; // the number of cartesian variables
    RefMolecularCoor mc_;
    RefMolecule mol_;

    RefSCVector cartesian_gradient_;
  protected:
    void failure(const char *);

    // this is just a wrapper around set_value()
    virtual void set_energy(double);

    // These are passed gradients and hessian in cartesian coordinates.
    // The _gradient and _hessian in internal coordinates are computed.
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    void x_to_molecule();
    void molecule_to_x();

  public:
    MolecularEnergy(const MolecularEnergy&);
    MolecularEnergy(const RefKeyVal&);
    MolecularEnergy(StateIn&);
    ~MolecularEnergy();

    void save_data_state(StateOut&);

    MolecularEnergy & operator=(const MolecularEnergy&);
    
    // a wrapper around value();
    virtual double energy();

    virtual RefMolecule molecule();
    virtual RefSCDimension moldim();
    
    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    void set_x(const RefSCVector&);

    RefSCVector get_cartesian_x();
    RefSCVector get_cartesian_gradient();

    virtual void print(ostream& = cout);
};
SavableState_REF_dec(MolecularEnergy);

class MolEnergyConvergence: public Convergence {
#   define CLASSNAME MolEnergyConvergence
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
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
    //will possibly grab the cartesian infomation if the Function
    //is a MolecularEnergy.
    void get_grad(const RefFunction &);
    void get_x(const RefFunction &);
    void get_nextx(const RefFunction &);

    // Return nonzero if the optimization has converged.
    int converged();
};
SavableState_REF_dec(MolEnergyConvergence);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
