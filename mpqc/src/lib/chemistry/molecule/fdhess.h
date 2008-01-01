//
// fdhess.h
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

#ifndef _chemistry_molecule_fdhess_h
#define _chemistry_molecule_fdhess_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream>

#include <chemistry/molecule/hess.h>
#include <chemistry/molecule/energy.h>

namespace sc {

/** Computes the molecular hessian by finite displacements of gradients.
    This will use the minimum number of displacements, each in the
    highest possible point group. */
class FinDispMolecularHessian: public MolecularHessian {
  protected:
    Ref<MolecularEnergy> mole_;
    // In case molecule must be given in lower symmetry, its actual
    // symmetry and the symmetry used to compute displacements is this
    Ref<PointGroup> displacement_point_group_;
    // The molecule's original point group for restoration at the end.
    Ref<PointGroup> original_point_group_;
    // The molecule's original geometry for restoration at the end and
    //computing displacements.
    RefSCVector original_geometry_;
    // the cartesian displacement size in bohr
    double disp_;
    // the accuracy for gradient calculations
    double accuracy_;
    // the number of completed displacements
    int ndisp_;
    // the number of irreps in the displacement point group
    int nirrep_;
    // whether or not to attempt a restart
    int restart_;
    // the name of the restart file
    std::string restart_file_;
    // whether or not to checkpoint
    int checkpoint_;
    // the name of the checkpoint file
    std::string checkpoint_file_;
    // only do the totally symmetric displacements
    int only_totally_symmetric_;
    // eliminate the cubic terms by doing an extra displacement for
    //each of the totally symmetry coordinates
    int eliminate_cubic_terms_;
    // use the gradient at the initial geometry to remove first order terms
    // (important if not at equilibrium geometry)
    int do_null_displacement_;
    // print flag
    int debug_;
    // a basis for the symmetrized cartesian coordinates
    RefSCMatrix symbasis_;
    // the gradients at each of the displacements
    RefSCVector *gradients_;

    void get_disp(int disp, int &irrep, int &index, double &coef);
    void do_hess_for_irrep(int irrep,
                           const RefSymmSCMatrix &dhessian,
                           const RefSymmSCMatrix &xhessian);
    void init();
    void restart();
  public:
    FinDispMolecularHessian(const Ref<MolecularEnergy>&);
    /** The FinDispMolecularHessian KeyVal constructor is used to generate a
        FinDispMolecularHessian object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>energy</tt><td>MolecularEnergy<td>none<td>This gives an
        object which will be used to compute the gradients needed to form
        the hessian.  If this is not specified, the object using
        FinDispMolecularHessian will, in some cases, fill it in
        appropriately.  However, even in these cases, it may be desirable
        to specify this keyword.  For example, this could be used in an
        optimization to compute frequencies using a lower level of theory.

        <tr><td><tt>debug</tt><td>boolean<td>false<td>If true,
        print out debugging information.

        <tr><td><tt>point_group</tt><td>PointGroup<td>none<td>
        The point group to use for generating the displacements.

        <tr><td><tt>restart</tt><td>boolean<td>true<td>If true, and a
        checkpoint file exists, restart from that file.

        <tr><td><tt>restart_file</tt><td>string
        <td><em>basename</em><tt>.ckpt.hess</tt><td>The name of
        the file where checkpoint information is written to or read from.

        <tr><td><tt>checkpoint</tt><td>boolean<td>true<td>If true,
        checkpoint intermediate data.

        <tr><td><tt>only_totally_symmetric</tt><td>boolean<td>false
        <td>If true, only follow totally symmetric displacments.  The
        hessian will not be complete, but it has enough information
        to use it in a geometry optimization.

        <tr><td><tt>eliminate_cubic_terms</tt><td>boolean<td>true<td>
        If true, then cubic terms will be eliminated.  This requires
        that two displacements are done for each totally symmetric
        coordinate, rather than one.  Setting this to false will reduce
        the accuracy, but the results will still probably be accurate
        enough for a geometry optimization.

        <tr><td><tt>do_null_displacement</tt><td>boolean<td>true<td>Run
        the calculation at the given geometry as well.

        <tr><td><tt>displacement</tt><td>double<td>1.0e-2<td>The size of
        the displacement in Bohr.

        <tr><td><tt>gradient_accuracy</tt><td>double<td><tt>displacement</tt>
        / 1000<td>The accuracy to which the gradients will be computed.

        </table>
    */
    FinDispMolecularHessian(const Ref<KeyVal>&);
    FinDispMolecularHessian(StateIn&);
    ~FinDispMolecularHessian();
    void save_data_state(StateOut&);

    /** These members are used to compute a cartesian hessian from
        gradients at finite displacements. */
    RefSymmSCMatrix compute_hessian_from_gradients();
    int ndisplace() const;
    int ndisplacements_done() const { return ndisp_; }
    RefSCMatrix displacements(int irrep) const;
    void displace(int disp);
    void original_geometry();
    void set_gradient(int disp, const RefSCVector &grad);
    void checkpoint_displacements(StateOut&);
    void restore_displacements(StateIn&);

    /** This returns the cartesian hessian.  If it has not yet been
        computed, it will be computed by finite displacements. */
    RefSymmSCMatrix cartesian_hessian();

    /// Set checkpoint option.
    void set_checkpoint(int c) { checkpoint_ = c; }
    /// Return the current value of the checkpoint option.
    int checkpoint() const { return checkpoint_; }

    void set_energy(const Ref<MolecularEnergy> &energy);
    MolecularEnergy* energy() const;

    Ref<SCMatrixKit> matrixkit() const { return mole_->matrixkit(); }
    RefSCDimension d3natom() const { return mole_->moldim(); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
