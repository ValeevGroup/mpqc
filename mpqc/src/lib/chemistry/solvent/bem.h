//
// bem.h
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

#ifndef _chemistry_solvent_bem_h
#define _chemistry_solvent_bem_h

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/isosurf/volume.h>
#include <math/isosurf/surf.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>

// This represents a solvent by a polarization charge on a dielectric
// boundary surface.
class BEMSolvent: public DescribedClass {
#   define CLASSNAME BEMSolvent
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefMolecule solute_;
    double dielectric_constant_;
    RefSCMatrixKit matrixkit_;
    RefSCMatrix system_matrix_i_;
    double f_;

    RefTriangulatedImplicitSurface surf_;

    double** alloc_array(int n, int m);
    void free_array(double**);

    // This holds the area associated with each vertex.  It is used
    // to convert charges to charge densities and back.
    double* vertex_area_;

    // Given charges compute surface charge density.
    void charges_to_surface_charge_density(double *charges);

    // Given surface charge density compute charges.
    void surface_charge_density_to_charges(double *charges);
  public:
    BEMSolvent(const RefKeyVal&);
    virtual ~BEMSolvent();

    // This should be called after everything is setup--the
    // molecule has the correct the geometry and all of the
    // parameters have been adjusted.
    void init();
    // This gets rid of the system matrix inverse.
    void done();

    int ncharge() { return surf_->nvertex(); }

    // NOTE: call allocation routines after init and free routines before done
    double** alloc_charge_positions() { return alloc_array(ncharge(), 3); }
    void free_charge_positions(double**a) { free_array(a); }

    double** alloc_normals()  { return alloc_array(ncharge(), 3); }
    void free_normals(double**a) { free_array(a); }

    double* alloc_efield_dot_normals()  { return new double[ncharge()]; }
    void free_efield_dot_normals(double*a) { delete[] a; }

    double* alloc_charges() { return new double[ncharge()]; }
    void free_charges(double*a) { delete[] a; }

    void charge_positions(double**);
    void normals(double**);

    // Given the efield dotted with the normals at the charge positions this
    // will compute a new set of charges.
    void compute_charges(double* efield_dot_normals, double* charge);

    // Given a set of charges and a total charge, this will normalize
    // the integrated charge to the charge that would be expected on
    // the surface if the given total charge were enclosed within it.
    void normalize_charge(double enclosed_charge, double* charges);

    // Given charges compute the interaction energy between the nuclei
    // and the point charges.
    double nuclear_interaction_energy(double** charge_positions,
                                      double* charge);

    // Given charges compute the interaction energy for just the surface.
    double self_interaction_energy(double** charge_positions, double *charge);
    
    // Given the charges, return the total polarization charge on the surface.
    double polarization_charge(double* charge);

    // this never needs to be called explicitly, but is here now for debugging
    void init_system_matrix();

    RefSCMatrixKit matrixkit() { return matrixkit_; }
};
DescribedClass_REF_dec(BEMSolvent);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
