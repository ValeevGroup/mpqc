//
// gaussbas.h
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

#ifndef _chemistry_qc_basis_gaussbas_h
#define _chemistry_qc_basis_gaussbas_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <util/state/state.h>
#include <util/state/array.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>

class GaussianShell;
class RefKeyVal;
class BasisFileSet;

SavableState_REF_fwddec(Integral)

class CartesianIter;
class SphericalTransformIter;

class GaussianBasisSet: public SavableState
{
#   define CLASSNAME GaussianBasisSet
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char* name_;
    GaussianShell** shell_;
    Arrayint shell_to_function_;
    Arrayint function_to_shell_;

    RefMolecule molecule_;

    RefSCMatrixKit matrixkit_;
    RefSCMatrixKit so_matrixkit_;
    RefSCDimension basisdim_;

    // these are needed for the routines to compute basis set values
    // they must be initialized with set_integral()
    CartesianIter **civec_;
    SphericalTransformIter **sivec_;

    int ncenter_;
    SSBArrayint shell_to_center_;
    SSBArrayint center_to_shell_;
    SSBArrayint center_to_nshell_;
    SSBArrayint center_to_nbasis_;

    int nshell_;
    int nbasis_;
    int nprim_;

    void recursively_get_shell(int&,RefKeyVal&,
                               const char*,const char*,BasisFileSet&,
                               int,int,int);

    void init(RefMolecule&,RefKeyVal&,
              BasisFileSet&,
              int have_userkeyval,
              int pure);
    void init2(int skip_ghosts=0);
    
  protected:
    GaussianBasisSet(const GaussianBasisSet&);
    virtual void set_matrixkit(const RefSCMatrixKit&);
    
  public:
    GaussianBasisSet(const RefKeyVal&);
    GaussianBasisSet(StateIn&);
    virtual ~GaussianBasisSet();

    void save_data_state(StateOut&);

    const char* name() const { return name_; }

    RefMolecule molecule() const { return molecule_; }
    RefSCMatrixKit matrixkit() { return matrixkit_; }
    RefSCMatrixKit so_matrixkit() { return so_matrixkit_; }
    RefSCDimension basisdim() { return basisdim_; }

    int ncenter() const;
    int nshell() const { return nshell_; }
    int nshell_on_center(int icenter) const;
    int shell_on_center(int icenter, int shell) const;
    int shell_to_center(int ishell) const { return shell_to_center_(ishell); }
    int nbasis() const { return nbasis_; }
    int nbasis_on_center(int icenter) const;
    int nprimitive() const { return nprim_; }

    int max_nfunction_in_shell() const;
    int max_ncartesian_in_shell(int aminc=0) const;
    int max_angular_momentum() const;
    // This is only need by integrals routines to set up
    // intermediate arrays.
    int max_ncontraction() const;
    int max_am_for_contraction(int con) const;
    int max_cartesian() const;

    int shell_to_function(int i) const { return shell_to_function_(i); }
    int function_to_shell(int i) const;

    // access to shells thru overall shell number
    const GaussianShell& operator()(int i) const { return *shell_[i]; }
    GaussianShell& operator()(int i) { return *shell_[i]; }
    const GaussianShell& operator[](int i) const { return *shell_[i]; }
    GaussianShell& operator[](int i) { return *shell_[i]; }
    const GaussianShell& shell(int i) const { return *shell_[i]; }
    GaussianShell& shell(int i) { return *shell_[i]; }

    // access to shells thru center number and relative shell number
    const GaussianShell& operator()(int icenter,int ishell) const;
    GaussianShell& operator()(int icenter,int ishell);
    const GaussianShell& shell(int i,int j) const { return operator()(i,j); }
    GaussianShell& shell(int i,int j) { return operator()(i,j); }

    // access to r thru center number
    double r(int icenter,int xyz) const;
    
    // compute the value for this basis set at position r
    // basis_values must be vector of length nbasis 
    int values(const SCVector3& r, double* basis_values) const;
    // g_values must be vector of length 3*nbasis
    // the data will be written in the order bf1_x, bf1_y, bf1_z, ...
    int grad_values(const SCVector3& r,
                    double*g_values,double* basis_values=0) const;
    // h_values must be vector of length 6*nbasis
    // the data will be written in the order bf1_xx, bf1_yx, bf1_yy,
    // bf1_zx, bf1_zy, bf1_zz, ...
    int hessian_values(const SCVector3& r, double *h_values,
                       double*g_values=0,double* basis_values=0) const;
    // this must be called before the above two routines to initialize
    // iterators that know the basis function order
    void set_integral(const RefIntegral&);

    // fill in matrix with a matrix that orthogonalizes the basis functions
    void ortho(const RefIntegral&, const RefSCMatrix&ortho);
    void ortho(const RefIntegral&,
               const RefSCMatrix&ortho, const RefSCMatrix&ortho_inverse);

    void print_brief(ostream& =cout) const;
    void print(ostream& =cout) const;
};

SavableState_REF_dec(GaussianBasisSet);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
