//
// deriv.h
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

#ifndef _chemistry_molecule_deriv_h
#define _chemistry_molecule_deriv_h

#include <iostream>

#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

namespace sc {

class MolecularEnergy;

/// @addtogroup ChemistryMolecule
/// @{

/** MolecularHessian is an abstract class that computes a molecule's second
    derivatives of the energy with respect to changes in the nuclear
    coordinates. */
class MolecularHessian: virtual public SavableState {
    static double desired_accuracy_default_;
    double desired_accuracy_;
    bool desired_accuracy_set_to_default_;
  protected:
    Ref<Molecule> mol_;
    RefSCDimension d3natom_;
    Ref<SCMatrixKit> matrixkit_;
  public:
    MolecularHessian();
    /** The MolecularHessian KeyVal constructor is used to generate a
        MolecularHessian derivative object from the input.  It reads the
        keywords below.

        <table border="1">
        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>molecule</tt><td>Molecule<td>none<td>The Molecule object.
        <tr><td><tt>accuracy</tt><td>real<td>1e-5<td>The desired accuracy of the hessian.
        </table>
    */
    MolecularHessian(const Ref<KeyVal>&);
    MolecularHessian(StateIn&);
    ~MolecularHessian();
    void save_data_state(StateOut&);

    RefSCDimension d3natom();
    Ref<SCMatrixKit> matrixkit() const { return matrixkit_; }

    /// Return the cartesian hessian.
    virtual RefSymmSCMatrix cartesian_hessian() = 0;

    /** Some MolecularHessian specializations require a molecular energy
        object.  The default implementations of this ignores the
        argument. */
    virtual void set_energy(const Ref<MolecularEnergy> &energy);
    /** This returns a MolecularEnergy object, if used by
        this specialization. Otherwise null is returned.  */
    virtual MolecularEnergy* energy() const;

    /** Find transformation matrix from cartesian to symmetry
        coordinates. */
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

    /**
     * Sets the desired accuracy
     * @param acc the desired accuracy
     */
    virtual void set_desired_accuracy(double acc);
    /**
     * Reports the desired accuracy
     * @return the desired accuracy
     */
    virtual double desired_accuracy() const;
    /**
     * @return whether the desired accuracy was set to default value
     */
    bool desired_accuracy_set_to_default() const { return desired_accuracy_set_to_default_; }
};


/** ReadMolecularHessian is an implementation of MolecularHessian
    that reads the hessian from a file. */
class ReadMolecularHessian: public MolecularHessian {
  protected:
    std::string filename_;
  public:
    /** The ReadMolecularHessian KeyVal constructor is used to generate a
        ReadMolecularHessian object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>filename</tt><td>string<td><em>basename</em>
        <tt>.hess</tt><td>The name of the file from which the hessian is
        read.

        </table>
    */
    ReadMolecularHessian(const Ref<KeyVal>&);
    ReadMolecularHessian(StateIn&);
    ~ReadMolecularHessian();
    void save_data_state(StateOut&);

    /// Return the hessian in cartesian coordinates.
    RefSymmSCMatrix cartesian_hessian();
};

/** GuessMolecularHessian is an implementation of MolecularHessian
    that estimates the hessian based on the internal coordinates. */
class GuessMolecularHessian: public MolecularHessian {
  protected:
    Ref<MolecularCoor> coor_;
  public:
    /** The GuessMolecularHessian KeyVal constructor is used to generate a
        GuessMolecularHessian object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>coor</tt><td>MolecularCoor<td>none<td>This gives
        the MolecularCoor object that is used to generate the guess
        hessian.  It does not have to be the same MolecularCoor
        object that is used to optimize the molecule.

        </table>
    */
    GuessMolecularHessian(const Ref<KeyVal>&);
    GuessMolecularHessian(StateIn&);
    ~GuessMolecularHessian();
    void save_data_state(StateOut&);

    /// Return the hessian in cartesian coordinates.
    RefSymmSCMatrix cartesian_hessian();
};

/** DiagMolecularHessian is an implementation of MolecularHessian
    that returns a hessian that is a diagonal matrix. */
class DiagMolecularHessian: public MolecularHessian {
  protected:
    double diag_;
  public:
    /** The DiagMolecularHessian KeyVal constructor is used to generate a
        DiagMolecularHessian object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>diag</tt><td>double<td>1.0<td>Specifies the diagonal
        elements of the hessian.

        </table>
    */
    DiagMolecularHessian(const Ref<KeyVal>&);
    DiagMolecularHessian(StateIn&);
    ~DiagMolecularHessian();
    void save_data_state(StateOut&);

    /// Return the hessian in cartesian coordinates.
    RefSymmSCMatrix cartesian_hessian();
};

////

/** MolecularGradient is an abstract class that computes a molecule's first
    derivatives of the energy with respect to changes in the nuclear
    coordinates. */
class MolecularGradient: virtual public SavableState {
    double desired_accuracy_;
  protected:
    Ref<Molecule> mol_;
    RefSCDimension d3natom_;
    Ref<SCMatrixKit> matrixkit_;
  public:
    MolecularGradient();
    /** The MolecularGradient KeyVal constructor is used to generate a
        MolecularGradient  object from the input.  It reads the
        keywords below.

        <table border="1">
        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>molecule</tt><td>Molecule<td>none<td>The Molecule object.
        <tr><td><tt>accuracy</tt><td>real<td>1e-4<td>The desired accuracy of the gradient.
        </table>
    */
    MolecularGradient(const Ref<KeyVal>&);
    MolecularGradient(StateIn&);
    ~MolecularGradient();
    void save_data_state(StateOut&);

    RefSCDimension d3natom();
    Ref<SCMatrixKit> matrixkit() const { return matrixkit_; }

    /// Return the cartesian hessian.
    virtual RefSCVector cartesian_gradient() = 0;

    /** Some MolecularGradient specializations require a molecular energy
        object.  The default implementations of this ignores the
        argument. */
    virtual void set_energy(const Ref<MolecularEnergy> &energy);
    /** This returns a MolecularEnergy object, if used by
        this specialization. Otherwise null is returned.  */
    virtual MolecularEnergy* energy() const;

    /// Write the gradient in a simple text format.
    static void write_cartesian_gradient(const char *filename,
                                         const Ref<Molecule> &m,
                                         const RefSCVector &grad);

    /// Read the hessian from a simple text format.
    static void read_cartesian_gradient(const char *filename,
                                        const Ref<Molecule> &m,
                                        const RefSCVector &grad);

    /**
     * Sets the desired accuracy
     * @param acc the desired accuracy
     */
    virtual void set_desired_accuracy(double acc);
    /**
     * Reports the desired accuracy
     * @return the desired accuracy
     */
    virtual double desired_accuracy() const;
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
