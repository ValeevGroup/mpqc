//
// gpetite.h --- definition of the generalized petite list class
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

#ifndef _chemistry_qc_basis_gpetite_h
#define _chemistry_qc_basis_gpetite_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdexcept>

#include <scconfig.h>
#include <util/misc/scint.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>

namespace sc {

/** If the shell loop structure has 8 fold symmetry, then
 *  this should be used as the template argument to GPetite4.
*/
class canonical_aaaa {
  public:
    canonical_aaaa();
    canonical_aaaa(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      long kl = (k>l?(((k*long(k+1))>>1)+l):(((l*long(l+1))>>1)+k));
      sc_int_least64_t
          off = (ij>kl?(((ij*sc_int_least64_t(ij+1))>>1)+kl)
                 :(((kl*sc_int_least64_t(kl+1))>>1)+ij));
      return off;
    }
};

/** If the shell loop structure has 2 fold symmetry between the first
 *  two indices, then this should be used as the template argument to
 *  GPetite4.
*/
class canonical_aabc {
    long nk_, nl_;
  public:
    canonical_aabc(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      return k + nk_*sc_int_least64_t(l + nl_*ij);
    }
};

/** If the shell loop structure has 2 fold symmetry between the last
 *  two indices, then this should be used as the template argument to
 *  GPetite4.
*/
class canonical_abcc {
    long nj_, nkl_;
  public:
    canonical_abcc(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long kl = (k>l?(((k*long(k+1))>>1)+l):(((l*long(l+1))>>1)+k));
      return kl + nkl_*sc_int_least64_t(j + nj_*i);
    }
};

/** If the shell loop structure has 2 fold symmetry between the first two
 *  indices and a 2 fold symmetry between the last two indices, then this
 *  should be used as the template argument to GPetite4.
*/
class canonical_aabb {
    long nij_;
  public:
    canonical_aabb(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      long kl = (k>l?(((k*long(k+1))>>1)+l):(((l*long(l+1))>>1)+k));
      return ij + nij_*sc_int_least64_t(kl);
    }
};

/** If the shell loop structure has 2 fold symmetry between the bra and the ket
 *  then this should be used as the template argument to GPetite4.
*/
class canonical_abab {
    long nj_;
  public:
    canonical_abab(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = i*nj_ + j;
      long kl = k*nj_ + l;
      sc_int_least64_t
          off = (ij>kl?(((ij*sc_int_least64_t(ij+1))>>1)+kl)
                 :(((kl*sc_int_least64_t(kl+1))>>1)+ij));
      return off;
    }
};

/** If the shell loop structure has no symmetry, then
 *  this should be used as the template argument to GPetite4.
*/
class canonical_abcd {
    int ni_, nj_, nk_;
  public:
    canonical_abcd(const Ref<GaussianBasisSet> bi,
                   const Ref<GaussianBasisSet> bj,
                   const Ref<GaussianBasisSet> bk,
                   const Ref<GaussianBasisSet> bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      return (i + ni_*sc_int_least64_t(j + nj_*long(k + nk_*l)));
    }
};

/** This class is an abstract base to a generalized four index petite list.
*/
class GenPetite4: public RefCount {
protected:
  bool c1_;
  int ng_;
  int **shell_map_i_;
  int **shell_map_j_;
  int **shell_map_k_;
  int **shell_map_l_;
  Ref<GaussianBasisSet> b1_, b2_, b3_, b4_;
public:
  GenPetite4(const Ref<GaussianBasisSet> &b1,
             const Ref<GaussianBasisSet> &b2,
             const Ref<GaussianBasisSet> &b3,
             const Ref<GaussianBasisSet> &b4);
  ~GenPetite4();
  virtual int in_p4(int i, int j, int k, int l) = 0;
};

/** This is a "factory" that prodces generalized four index petite list objects.
*/
extern Ref<GenPetite4>
construct_gpetite(const Ref<GaussianBasisSet> &b1,
                  const Ref<GaussianBasisSet> &b2,
                  const Ref<GaussianBasisSet> &b3,
                  const Ref<GaussianBasisSet> &b4);

/** This class provides a generalized four index petite list.
 * The template argument is a class that computes an canonical offset
 * given four indices for the particular shell loop structure employed.
 * Example template class parameters are canonical_aaaa, canonical_aabc,
 * canonical_aabb, and canonical_abcd.
 */
template <class C4>
class GPetite4: public GenPetite4 {
    C4 c_;
  public:
    GPetite4(const Ref<GaussianBasisSet> &b1,
             const Ref<GaussianBasisSet> &b2,
             const Ref<GaussianBasisSet> &b3,
             const Ref<GaussianBasisSet> &b4,
             const C4& c);
    ~GPetite4();
    int in_p4(int i, int j, int k, int l);
};

template <class C4>
inline int
GPetite4<C4>::in_p4(int i, int j, int k, int l)
{
  if (c1_) return 1;

  sc_int_least64_t ijkl = c_.offset(i,j,k,l);
  int nijkl = 1;

  for (int g=1; g < ng_; g++) {
      int gi = shell_map_i_[i][g];
      int gj = shell_map_j_[j][g];
      int gk = shell_map_k_[k][g];
      int gl = shell_map_l_[l][g];
      sc_int_least64_t gijkl = c_.offset(gi,gj,gk,gl);

      if (gijkl > ijkl) return 0;
      else if (gijkl == ijkl) nijkl++;
  }

  return ng_/nijkl;
}

}



#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
