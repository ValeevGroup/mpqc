//
// hess.h
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

#ifndef _chemistry_molecule_hess_h
#define _chemistry_molecule_hess_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream>

#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

namespace sc {

class MolecularEnergy;

class MolecularHessian: virtual public SavableState {
  protected:
    Ref<Molecule> mol_;
    RefSCDimension d3natom_;
    Ref<SCMatrixKit> matrixkit_;
  public:
    MolecularHessian();
    MolecularHessian(const Ref<KeyVal>&);
    MolecularHessian(StateIn&);
    ~MolecularHessian();
    void save_data_state(StateOut&);

    RefSCDimension d3natom();
    Ref<SCMatrixKit> matrixkit() const { return matrixkit_; }

    // Return the cartesian hessian.
    virtual RefSymmSCMatrix cartesian_hessian() = 0;

    // Some MolecularHessian specializations require a molecular
    //energy object.  The default implementations of these ignore
    //the argument or return null.
    virtual void set_energy(const Ref<MolecularEnergy> &energy);
    virtual MolecularEnergy* energy() const;

    // Find transformation matrix from cartesian to symmetry
    // coordinates.
    static RefSCMatrix cartesian_to_symmetry(const Ref<Molecule> &m,
                                             Ref<PointGroup> pg = 0,
                                             Ref<SCMatrixKit> kit = 0);

    /// Write the hessian in a simple text format.
    static void write_cartesian_hessian(const char *filename,
                                        const Ref<Molecule> &m,
                                        const RefSymmSCMatrix &hess);

    /// Read the hessian from a simple text format.
    static void read_cartesian_hessian(const char *filename,
                                       const Ref<Molecule> &m,
                                       const RefSymmSCMatrix &hess);
};


class ReadMolecularHessian: public MolecularHessian {
  protected:
    char *filename_;
  public:
    ReadMolecularHessian(const Ref<KeyVal>&);
    ReadMolecularHessian(StateIn&);
    ~ReadMolecularHessian();
    void save_data_state(StateOut&);

    RefSymmSCMatrix cartesian_hessian();
};

class GuessMolecularHessian: public MolecularHessian {
  protected:
    Ref<MolecularCoor> coor_;
  public:
    GuessMolecularHessian(const Ref<KeyVal>&);
    GuessMolecularHessian(StateIn&);
    ~GuessMolecularHessian();
    void save_data_state(StateOut&);

    RefSymmSCMatrix cartesian_hessian();
};

class DiagMolecularHessian: public MolecularHessian {
  protected:
    double diag_;
  public:
    DiagMolecularHessian(const Ref<KeyVal>&);
    DiagMolecularHessian(StateIn&);
    ~DiagMolecularHessian();
    void save_data_state(StateOut&);

    RefSymmSCMatrix cartesian_hessian();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
