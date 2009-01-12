//
// petite.h --- definition of the PetiteList class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _chemistry_qc_basis_petite_h
#define _chemistry_qc_basis_petite_h

#ifdef __GNUC__
#pragma interface
#endif

#include <scconfig.h>
#include <iostream>
#include <scconfig.h>

#include <util/misc/scint.h>
#include <util/ref/ref.h>
#include <math/scmat/blocked.h>
#include <math/scmat/offset.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/integral.h>

// //////////////////////////////////////////////////////////////////////////

namespace sc {

inline sc_int_least64_t
ij_offset64(sc_int_least64_t i, sc_int_least64_t j)
{
  return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
}

inline sc_int_least64_t
i_offset64(sc_int_least64_t i)
{
  return ((i*(i+1)) >> 1);
}

// //////////////////////////////////////////////////////////////////////////
// These are helper functions for PetiteList and GenericPetiteList4

int **compute_atom_map(const Ref<GaussianBasisSet> &);
void delete_atom_map(int **atom_map, const Ref<GaussianBasisSet> &);

int **compute_shell_map(int **atom_map, const Ref<GaussianBasisSet> &);
void delete_shell_map(int **shell_map, const Ref<GaussianBasisSet> &);

// //////////////////////////////////////////////////////////////////////////

struct contribution {
    int bfn;
    double coef;

    contribution();
    contribution(int b, double c);
    ~contribution();
};

struct SO {
    int len;
    int length;
    contribution *cont;

    SO();
    SO(int);
    ~SO();

    SO& operator=(const SO&);

    void set_length(int);
    void reset_length(int);

    // is this equal to so to within a sign
    int equiv(const SO& so);
};

struct SO_block {
    int len;
    SO *so;

    SO_block();
    SO_block(int);
    ~SO_block();

    void set_length(int);
    void reset_length(int);

    int add(SO& s, int i);
    void print(const char *title);
};

// //////////////////////////////////////////////////////////////////////////

class PetiteList : public RefCount {
  private:
    int natom_;
    int nshell_;
    int ng_;
    int nirrep_;
    int nblocks_;
    int c1_;

    Ref<GaussianBasisSet> gbs_;
    Ref<Integral> ints_;

    char *p1_;        // p1[n] is 1 if shell n is in the group P1
    int **atom_map_;  // atom_map[n][g] is the atom that symop g maps atom n
                      // into
    int **shell_map_; // shell_map[n][g] is the shell that symop g maps shell n
                      // into
    char *lamij_;     // see Dupuis & King, IJQC 11,613,(1977)

    int *nbf_in_ir_;

    void init();

  public:
    PetiteList(const Ref<GaussianBasisSet>&, const Ref<Integral>&);
    ~PetiteList();

    Ref<GaussianBasisSet> basis() { return gbs_; }
    Ref<Integral> integral() { return ints_; }
    Ref<PetiteList> clone() { return new PetiteList(gbs_, ints_); }

    int nirrep() const { return nirrep_; }
    int order() const { return ng_; }
    int atom_map(int n, int g) const { return (c1_) ? n : atom_map_[n][g]; }
    int shell_map(int n, int g) const { return (c1_) ? n : shell_map_[n][g]; }
    int lambda(int ij) const { return (c1_) ? 1 : (int) lamij_[ij]; }
    int lambda(int i, int j) const
                          { return (c1_) ? 1 : (int) lamij_[ij_offset(i,j)]; }

    int in_p1(int n) const { return (c1_) ? 1 : (int) p1_[n]; }
    int in_p2(int ij) const { return (c1_) ? 1 : (int) lamij_[ij]; }
    /// Same as previous, except for it takes i and j separately
    int in_p2(int i, int j) const { return (c1_) ? 1 : (int) lamij_[ij_offset(i,j)]; }
    int in_p4(int ij, int kl, int i, int j, int k, int l) const;
    /// Same as previous, except for doesn't assume ij > kl and recomputes them
    int in_p4(int i, int j, int k, int l) const;

    int nfunction(int i) const
                            { return (c1_) ? gbs_->nbasis() : nbf_in_ir_[i]; }

    int nblocks() const { return nblocks_; }

    void print(std::ostream& =ExEnv::out0(), int verbose=1);

    // these return blocked dimensions
    RefSCDimension AO_basisdim();
    RefSCDimension SO_basisdim();

    // return the basis function rotation matrix R(g)
    RefSCMatrix r(int g);

    // return information about the transformation from AOs to SOs
    SO_block * aotoso_info();
    RefSCMatrix aotoso();
    RefSCMatrix sotoao();

    // given a skeleton matrix, form the symmetrized matrix in the SO basis
    void symmetrize(const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym);

    // transform a matrix from AO->SO or SO->AO.
    // this can take either a blocked or non-blocked AO basis matrix.
    RefSymmSCMatrix to_SO_basis(const RefSymmSCMatrix&);

    // this returns a non-blocked AO basis matrix.
    RefSymmSCMatrix to_AO_basis(const RefSymmSCMatrix&);

    // these two are just for eigenvectors
    // returns non-blocked AO basis eigenvectors
    RefSCMatrix evecs_to_AO_basis(const RefSCMatrix&);
    // returns blocked SO basis eigenvectors
    RefSCMatrix evecs_to_SO_basis(const RefSCMatrix&);
};

inline int
PetiteList::in_p4(int ij, int kl, int i, int j, int k, int l) const
{
  if (c1_)
    return 1;

  sc_int_least64_t ijkl = i_offset64(ij)+kl;
  int nijkl=1;

  for (int g=1; g < ng_; g++) {
    int gij = ij_offset(shell_map_[i][g],shell_map_[j][g]);
    int gkl = ij_offset(shell_map_[k][g],shell_map_[l][g]);
    sc_int_least64_t gijkl = ij_offset64(gij,gkl);

    if (gijkl > ijkl)
      return 0;
    else if (gijkl == ijkl)
      nijkl++;
  }

  return ng_/nijkl;
}

inline int
PetiteList::in_p4(int i, int j, int k, int l) const
{
  if (c1_)
    return 1;

  int ij = ij_offset(i,j);
  int kl = ij_offset(k,l);
  sc_int_least64_t ijkl = ij_offset64(ij,kl);
  int nijkl=1;

  for (int g=1; g < ng_; g++) {
    int gij = ij_offset(shell_map_[i][g],shell_map_[j][g]);
    int gkl = ij_offset(shell_map_[k][g],shell_map_[l][g]);
    sc_int_least64_t gijkl = ij_offset64(gij,gkl);

    if (gijkl > ijkl)
      return 0;
    else if (gijkl == ijkl)
      nijkl++;
  }

  return ng_/nijkl;
}

}



#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
