//
// molfreq.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_molfreq_h
#define _chemistry_qc_molfreq_h

#include <iostream.h>
#include <math/scmat/matrix.h>
#include <util/render/animate.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molrender.h>

DescribedClass_REF_fwddec(RenderedObject);
DescribedClass_REF_fwddec(RenderedMolecule);
DescribedClass_REF_fwddec(MolFreqAnimate);

class MolecularFrequencies: public SavableState {
#   define CLASSNAME MolecularFrequencies
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefMolecularEnergy mole_;
    RefSCMatrixKit kit_;
    RefSCMatrixKit symkit_;
    RefMolecule mol_;
    // In case molecule must be given in lower symmetry, its actual
    // symmetry and the symmetry used to compute displacements is this
    RefPointGroup displacement_point_group_;
    RefPointGroup original_point_group_;
    RefSCVector original_geometry_;
    // the cartesian displacement size in bohr
    double disp_;
    // displacements for each irrep
    int ndisp_;
    int nirrep_;
    RefSCMatrix *displacements_;
    RefSCVector *gradients_;

    // the number of external degrees of freedom
    int nexternal_;

    // the number of frequencies per irrep
    int *nfreq_;
    // the frequencies for each irrep
    double **freq_;
    // the symmetry blocked dimension for internal motions
    RefSCDimension disym_;
    // the normal coordinates
    RefSCMatrix normco_;
    // the cartesian dimension
    RefSCDimension d3natom_;
    // the blocked cartesian dimension
    RefSCDimension bd3natom_;

    void get_disp(int disp, int &irrep, int &index, double &coef);
    void do_freq_for_irrep(int irrep,
                           const RefDiagSCMatrix &m,
                           const RefSymmSCMatrix &dhessian,
                           const RefSymmSCMatrix &xhessian);
    int debug_;
  public:
    MolecularFrequencies(const RefKeyVal &);
    MolecularFrequencies(StateIn &);
    ~MolecularFrequencies();
    void save_data_state(StateOut&);

    void compute_displacements();
    void compute_frequencies_from_gradients();
    int ndisplace() const;
    void displace(int disp);
    void original_geometry();
    void set_gradient(int disp, const RefSCVector &grad);

    void thermochemistry(int degeneracy, double temp=298.15, double pres=1.0);

    void animate(const RefRender&, const RefMolFreqAnimate&);

    RefSCMatrix normal_coordinates() { return normco_; }

    RefSCMatrixKit matrixkit() { return kit_; }
    RefSCMatrixKit symmatrixkit() { return symkit_; }
    RefSCDimension d3natom() { return d3natom_; }
};

SavableState_REF_dec(MolecularFrequencies);

class MolFreqAnimate: public AnimatedObject {
#   define CLASSNAME MolFreqAnimate
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    RefRenderedMolecule renmol_;
    RefMolecularFrequencies molfreq_;
    RefMolecularEnergy dependent_mole_;
    int irrep_;
    int mode_;
    int nframe_;
  public:
    MolFreqAnimate(const RefKeyVal &);
    virtual ~MolFreqAnimate();

    void set_mode(int i, int j) { irrep_ = i; mode_ = j; }
    int nobject();
    RefRenderedObject object(int iobject);
};
DescribedClass_REF_dec(MolFreqAnimate);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
