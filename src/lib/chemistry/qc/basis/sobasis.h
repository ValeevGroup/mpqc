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

#include <chemistry/qc/basis/basis.h>

namespace sc {

/** SOTransformShell describes how an AO function contributes to an SO
    function in a particular SO shell. */
class SOTransformFunction {
  public:
    /// The coefficient of the AO.
    double coef;
    /// The AO function number.
    int aofunc;
    /// The SO function number.
    int sofunc;
    /// The SO function's irrep.
    int irrep;
};

/** SOTransformShell maintains a list of AO functions contribute to an SO
    function in a particular SO shell.  The information is stored in
    objects of type SOTransformFunction. */
class SOTransformShell {
  public:
    /// The number of the AO shell from which these functions come.
    int aoshell;
    /// The number of AO/SO function pairs contributing.
    int nfunc;
    /// The array of SOTransformFunction objects describing the transform.
    SOTransformFunction *func;
    SOTransformShell();
    ~SOTransformShell();
    /// Add another function to the transform.
    void add_func(int irrep, double coef, int aofunc, int sofunc);
};

/** SOTransform maintains a list of AO shells that are be used
    to compute the SO.  The information is stored in objects of
    type SOTransformShell. */
class SOTransform {
  public:
    int naoshell_allocated;
    /// The number of AO shells that make up this SO shell.
    int naoshell;
    /// The SOTransformShell object for each AO.
    SOTransformShell *aoshell;
    SOTransform();
    ~SOTransform();
    void set_naoshell(int n);
    /// Adds another term to the transform.
    void add_transform(int aoshell, int irrep,
                       double coef, int aofunc, int sofunc);
};

/** A SOBasis object describes the transformation from an atomic orbital
    basis to a symmetry orbital basis. The cental concept here is "SO shell".
    An SO shell is defined as a collection of symmetry-adapated basis functions that
    arise from one set of AO shells related by symmetry (an AO shell orbit).
    Thus the number of SO shells = the number of symmetry-unique AO shells,
    and each SO shell contains basis functions of potentially several symmetries.
    */
class SOBasis : public RefCount {
  protected:
    Ref<GaussianBasisSet> basis_;
    int nshell_;       //< number of unique AO shells
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
    /// Create an SOBasis object given a GaussianBasisSet and Integral objects.
    SOBasis(const Ref<GaussianBasisSet> &, const Ref<Integral>&);
    ~SOBasis();

    /// Return the number of shells.
    int nshell() const { return nshell_; }
    /// Return the number of irreps.
    int nirrep() const { return nirrep_; }
    int ncomponent(int iirrep) const { return ncomp_[iirrep]; }
    /// Return the number of functions in the given irrep.
    int nfunction_in_irrep(int irrep) const { return nfunc_in_irrep_[irrep]; }
    /// Return the offset for the first function of the given irrep.
    int function_offset_for_irrep(int irrep) const;
    /// Return the number of functions in the given shell.
    int nfunction(int ishell) const;
    /** Return the number of functions in the AO shell that make up
        the given SO shell. */
    int naofunction(int ishell) const { return naofunc_[ishell]; }
    /// Returns the number of functions in the shell in a given irrep.
    int nfunction(int ishell, int iirrep) const;
    /** Returns the maximum number of functions in a shell (summed over all
        irreps) */
    int max_nfunction_in_shell() const;
    /** Normally, SO shell numbering starts at zero within each irrep.
        This returns an offset to make SO shell numbers unique within the
        shell. */
    int function_offset_within_shell(int ishell, int iirrep) const;

    /** Convert the SO shell number to the overall number of the first
        function within that shell. */
    int function(int ishell);

    /// Convert SO shell and function number within shell to irrep.
    int irrep(int ishell, int ifunc) const;
    /// Convert SO shell and function number to number within irrep.
    int function_within_irrep(int ishell, int ifunc) const;

    /// Return the SOTransform object for the given shell.
    const SOTransform &trans(int ishell) const { return trans_[ishell]; }

    void print(std::ostream &o=ExEnv::out0()) const;
};


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
SOBasis::function_offset_for_irrep(int irrep) const
{
  int r = 0;
  for (int i=0; i<irrep; i++) {
      r += nfunc_in_irrep_[i];
    }
  return r;
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

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
