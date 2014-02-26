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

#ifndef _chemistry_qc_basis_gpetite_h
#define _chemistry_qc_basis_gpetite_h

#include <stdexcept>

#include <mpqc_config.h>
#include <util/misc/scint.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>

namespace sc {

/** If the shell loop structure has 8 fold symmetry, then
 *  this should be used as the template argument to GenericPetiteList4.
*/
class canonical_aaaa {
  public:
    canonical_aaaa();
    canonical_aaaa(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
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
 *  GenericPetiteList4.
*/
class canonical_aabc {
    long nk_, nl_;
  public:
    canonical_aabc(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      return k + nk_*sc_int_least64_t(l + nl_*ij);
    }
};

/** If the shell loop structure has 2 fold symmetry between the last
 *  two indices, then this should be used as the template argument to
 *  GenericPetiteList4.
*/
class canonical_abcc {
    long nj_, nkl_;
  public:
    canonical_abcc(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long kl = (k>l?(((k*long(k+1))>>1)+l):(((l*long(l+1))>>1)+k));
      return kl + nkl_*sc_int_least64_t(j + nj_*i);
    }
};

/** If the shell loop structure has 2 fold symmetry between the first two
 *  indices and a 2 fold symmetry between the last two indices, then this
 *  should be used as the template argument to GenericPetiteList4.
*/
class canonical_aabb {
    long nij_;
  public:
    canonical_aabb(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      long kl = (k>l?(((k*long(k+1))>>1)+l):(((l*long(l+1))>>1)+k));
      return ij + nij_*sc_int_least64_t(kl);
    }
};

/** If the shell loop structure has 2 fold symmetry between the bra and the ket
 *  then this should be used as the template argument to GenericPetiteList4.
*/
class canonical_abab {
    long nj_;
  public:
    canonical_abab(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
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
 *  this should be used as the template argument to GenericPetiteList4.
*/
class canonical_abcd {
    int ni_, nj_, nk_;
  public:
    canonical_abcd(const Ref<GaussianBasisSet>& bi,
                   const Ref<GaussianBasisSet>& bj,
                   const Ref<GaussianBasisSet>& bk,
                   const Ref<GaussianBasisSet>& bl
                   );
    sc_int_least64_t offset(int i, int j, int k, int l) {
      return (i + ni_*sc_int_least64_t(j + nj_*long(k + nk_*l)));
    }
};

/** This class is an abstract base to a generalized four index petite list.
*/
class GPetiteList4: public RefCount {
  bool c1_;
  int ng_;
protected:
  int **shell_map_i_;
  int **shell_map_j_;
  int **shell_map_k_;
  int **shell_map_l_;
  Ref<GaussianBasisSet> b1_, b2_, b3_, b4_;
  bool c1() const { return c1_; }
  int order() const { return ng_; }
public:
  GPetiteList4(const Ref<GaussianBasisSet> &b1,
             const Ref<GaussianBasisSet> &b2,
             const Ref<GaussianBasisSet> &b3,
             const Ref<GaussianBasisSet> &b4);
  ~GPetiteList4();
  const Ref<GaussianBasisSet>& basis1() const { return b1_; }
  const Ref<GaussianBasisSet>& basis2() const { return b2_; }
  const Ref<GaussianBasisSet>& basis3() const { return b3_; }
  const Ref<GaussianBasisSet>& basis4() const { return b4_; }
  const Ref<PointGroup>& point_group() const { return b1_->molecule()->point_group(); }
  virtual int in(int i, int j, int k, int l) = 0;
};

/** This class provides a generalized four index petite list.
 * The template argument is a class that computes an canonical offset
 * given four indices for the particular shell loop structure employed.
 * Example template class parameters are canonical_aaaa, canonical_aabc,
 * canonical_aabb, and canonical_abcd.
 */
template <class C4>
class GenericPetiteList4 : public GPetiteList4 {
    C4 c_;
  public:
    GenericPetiteList4(const Ref<GaussianBasisSet> &b1,
             const Ref<GaussianBasisSet> &b2,
             const Ref<GaussianBasisSet> &b3,
             const Ref<GaussianBasisSet> &b4) : GPetiteList4(b1,b2,b3,b4),
             c_(b1,b2,b3,b4) {}
    ~GenericPetiteList4() {}
    int in(int i, int j, int k, int l);
};

template <class C4>
inline int
GenericPetiteList4<C4>::in(int i, int j, int k, int l)
{
  if (c1()) return 1;

  sc_int_least64_t ijkl = c_.offset(i,j,k,l);
  int nijkl = 1;

  const int ng = order();
  for (int g=1; g < ng; g++) {
      int gi = shell_map_i_[i][g];
      int gj = shell_map_j_[j][g];
      int gk = shell_map_k_[k][g];
      int gl = shell_map_l_[l][g];
      sc_int_least64_t gijkl = c_.offset(gi,gj,gk,gl);

      if (gijkl > ijkl) return 0;
      else if (gijkl == ijkl) nijkl++;
  }

  return ng/nijkl;
}

/// Can be used as a template argument to GenericPetiteList2
class canonical_aa {
  public:
    canonical_aa(const Ref<GaussianBasisSet>& bi,
                 const Ref<GaussianBasisSet>& bj
    );
    sc_int_least64_t offset(int i, int j) {
      const long ij = (i>j?(((i*long(i+1))>>1)+j):(((j*long(j+1))>>1)+i));
      return ij;
    }
};

/** Can be used as a template argument to GenericPetiteList2
*/
class canonical_ab {
    int ni_;
  public:
    canonical_ab(const Ref<GaussianBasisSet>& bi,
                 const Ref<GaussianBasisSet>& bj
    );
    sc_int_least64_t offset(int i, int j) {
      return (i + ni_*sc_int_least64_t(j));
    }
};

/** This class is an abstract base to a generalized 2-index petite list.
*/
class GPetiteList2: public RefCount {
  bool c1_;
  int ng_;
protected:
  int **shell_map_i_;
  int **shell_map_j_;
  Ref<GaussianBasisSet> b1_, b2_;
  bool c1() const { return c1_; }
  int order() const { return ng_; }
public:
  GPetiteList2(const Ref<GaussianBasisSet> &b1,
               const Ref<GaussianBasisSet> &b2);
  ~GPetiteList2();
  const Ref<GaussianBasisSet>& basis1() const { return b1_; }
  const Ref<GaussianBasisSet>& basis2() const { return b2_; }
  const Ref<PointGroup>& point_group() const { return b1_->molecule()->point_group(); }
  virtual int in(int i, int j) = 0;
};

/** This class provides a generalized 2-index petite list.
 * The template argument is a class that computes an canonical offset
 * given two indices for the particular shell loop structure employed.
 * Example template class parameters are canonical_aa and canonical_ab.
 */
template <class C2>
class GenericPetiteList2 : public GPetiteList2 {
    C2 c_;
  public:
    GenericPetiteList2(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2) : GPetiteList2(b1,b2),
                       c_(b1,b2) {}
    ~GenericPetiteList2() {}
    int in(int i, int j);
    void symmetrize(const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym) const;
    void symmetrize(const RefSCMatrix& skel, const RefSCMatrix& sym) const;
};

template <class C2>
inline int
GenericPetiteList2<C2>::in(int i, int j)
{
  if (c1()) return 1;

  sc_int_least64_t ij = c_.offset(i,j);
  int nij = 1;

  const int ng = order();
  for (int g=1; g < ng; g++) {
      int gi = shell_map_i_[i][g];
      int gj = shell_map_j_[j][g];
      const sc_int_least64_t gij = c_.offset(gi,gj);

      if (gij > ij) return 0;
      else if (gij == ij) nij++;
  }

  return ng/nij;
}

/** Produces generalized 2 and 4-index petite list objects.
*/
struct GPetiteListFactory {
  static Ref<GPetiteList2> plist2(const Ref<GaussianBasisSet> &b1,
                                  const Ref<GaussianBasisSet> &b2);
  static Ref<GPetiteList4> plist4(const Ref<GaussianBasisSet> &b1,
                                  const Ref<GaussianBasisSet> &b2,
                                  const Ref<GaussianBasisSet> &b3,
                                  const Ref<GaussianBasisSet> &b4);
};

/// Uses plist2 to convert the "skeleton" matrix into the full matrix. Only applicable when the two basis sets are equivalent.
void symmetrize(const Ref<GPetiteList2>& plist2, const Ref<Integral>& integral, const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym);
/// Uses plist2 to convert the "skeleton" matrix into the full matrix.
void symmetrize(const Ref<GPetiteList2>& plist2, const Ref<Integral>& integral, const RefSCMatrix& skel, const RefSCMatrix& sym);
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
