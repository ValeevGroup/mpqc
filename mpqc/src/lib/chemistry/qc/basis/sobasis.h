//
// sobasis.h --- definition of the Integral class
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_basis_sobasis_h
#define _chemistry_qc_basis_sobasis_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/basis/basis.h>

class SOTransformFunction {
  public:
    double coef;
    int aofunc;
    int sofunc;
    int irrep;
};

class SOTransformShell {
  public:
    int aoshell;
    int nfunc;
    SOTransformFunction *func;
    SOTransformShell();
    ~SOTransformShell();
    void add_func(int irrep, double coef, int aofunc, int sofunc);
};

class SOTransform {
  public:
    int naoshell_allocated;
    int naoshell;
    SOTransformShell *aoshell;
    SOTransform();
    ~SOTransform();
    void set_naoshell(int n);
    void add_transform(int aoshell, int irrep,
                       double coef, int aofunc, int sofunc);
};

class SOBasis : public VRefCount {
  protected:
    RefGaussianBasisSet basis_;
    int nshell_;
    int nirrep_;
    int *ncomp_;
    int **nfunc_;
    int *naofunc_;
    int **funcoff_;

    int *nfunc_in_irrep_;
    int *func_;
    int *irrep_;
    int *func_within_irrep_;

    SOTransform *trans_;

  public:
    SOBasis(const RefGaussianBasisSet &, const RefIntegral&);
    ~SOBasis();

    int nshell() const { return nshell_; }
    int nirrep() const { return nirrep_; }
    int ncomponent(int iirrep) const { return ncomp_[iirrep]; }
    // returns the number of functions in the shell
    int nfunction(int ishell) const;
    // return the number of functions in the ao shell that make up
    // the given so shell
    int naofunction(int ishell) const { return naofunc_[ishell]; }
    // returns the number of functions in the shell in a given irrep
    int nfunction(int ishell, int iirrep) const;
    // returns the max number of functions in a shell (summed over all irreps)
    int max_nfunction_in_shell() const;
    // normally, so shell numbering starts at zero within each irrep
    // this returns an offset to make so shell numbers unique within the shell
    int function_offset_within_shell(int ishell, int iirrep) const;

    // Convert the so shell number to the overall number of the first
    // function within that shell.
    int function(int ishell);

    // so shell and function number with shell to irrep
    int irrep(int ishell, int ifunc) const;
    // so shell and function number to number within irrep
    int function_within_irrep(int ishell, int ifunc) const;

    const SOTransform &trans(int i) const { return trans_[i]; }

    void print(ostream &o=cout) const;
};
REF_dec(SOBasis);

inline int
SOBasis::function(int ishell)
{
  return func_[ishell];
}

inline int
SOBasis::irrep(int ishell, int ifunc) const
{
  return irrep_[func_[ishell]+ifunc];
}

inline int
SOBasis::function_within_irrep(int ishell, int ifunc) const
{
  return func_within_irrep_[func_[ishell]+ifunc];
}

inline int
SOBasis::nfunction(int ishell, int iirrep) const
{
  return nfunc_[ishell][iirrep];
}

inline int
SOBasis::function_offset_within_shell(int ishell, int iirrep) const
{
  return funcoff_[ishell][iirrep];
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
