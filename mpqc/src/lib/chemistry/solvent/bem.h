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

namespace sc {

// This represents a solvent by a polarization charge on a dielectric
// boundary surface.
class BEMSolvent: public DescribedClass {
  private:
    int debug_;

    Ref<Molecule> solute_;
    Ref<Molecule> solvent_;
    double solvent_density_;
    double dielectric_constant_;
    Ref<SCMatrixKit> matrixkit_;
    RefSCMatrix system_matrix_i_;
    double f_;
    Ref<MessageGrp> grp_;

    double area_;
    double volume_;
    double computed_enclosed_charge_;
    double edisp_;
    double erep_;

    Ref<TriangulatedImplicitSurface> surf_;

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
    BEMSolvent(const Ref<KeyVal>&);
    virtual ~BEMSolvent();

    // This should be called after everything is setup--the
    // molecule has the correct the geometry and all of the
    // parameters have been adjusted.
    void init();
    // This gets rid of the system matrix inverse and, optionally, the surface.
    void done(int clear_surface = 1);

    int ncharge() { return surf_->nvertex(); }

    Ref<Molecule> solvent() { return solvent_ ;}
    double solvent_density() { return solvent_density_ ;}

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

    // Given charges and nuclear charges compute their interation energy.
    double nuclear_charge_interaction_energy(double *nuclear_charge,
                                             double** charge_positions,
                                             double* charge);

    // Given charges compute the interaction energy between the nuclei
    // and the point charges.
    double nuclear_interaction_energy(double** charge_positions,
                                      double* charge);

    // Given charges compute the interaction energy for just the surface.
    double self_interaction_energy(double** charge_positions, double *charge);
    
    // Given the charges, return the total polarization charge on the surface.
    double polarization_charge(double* charge);

    // Return the area (available after compute_charges called).
    double area() const { return area_; }
    // Return the volume (available after compute_charges called).
    double volume() const { return volume_; }
    // Return the enclosed charge (available after compute_charges called).
    double computed_enclosed_charge() const {
      return computed_enclosed_charge_;
    }

    double disp() {return edisp_;}
    double rep()  {return erep_;}
    double disprep();

    // this never needs to be called explicitly, but is here now for debugging
    void init_system_matrix();

    Ref<TriangulatedImplicitSurface> surface() const { return surf_; }

    Ref<SCMatrixKit> matrixkit() { return matrixkit_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
