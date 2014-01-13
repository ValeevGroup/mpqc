//

// cadf_iters.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jan 10, 2014
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

#ifndef _chemistry_qc_scf_cadf_iters_h
#define _chemistry_qc_scf_cadf_iters_h

#include <chemistry/qc/scf/clhf.h>
#include <util/misc/property.h>
#include <utility>

#define DEFAULT_TARGET_BLOCK_SIZE 1000  // functions

namespace sc {

#define DUMP(expr) std::cout << #expr << " = " << (expr) << std::endl;
#define out_assert(a, op, b) assert(a op b || ((std::cout << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))

#define NOT_ASSIGNED -1

// Forward declarations
class ShellData;
class BasisFunctionData;
class ShellBlockData;

typedef enum {
  NoRestrictions = 0,
  SameCenter = 1,
  SameAngularMomentum = 2
} BlockCompositionRequirement;

//############################################################################//

namespace detail {

template <typename DataContainer>
class basis_iterator {

  protected:

    GaussianBasisSet* basis_;
    GaussianBasisSet* dfbasis_;
    int first_index_;
    int last_index_;

  public:

    enum { NoLastIndex = -1 };

    typedef DataContainer value_type;

    basis_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis),
        dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_index_(0),
        last_index_(-1)
    { }

    basis_iterator(
        const Ref<GaussianBasisSet>& basis,
        int last_index,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis),
        dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_index_(0),
        last_index_(last_index)
    { }

    basis_iterator(
        const Ref<GaussianBasisSet>& basis,
        int first_index,
        int last_index,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis),
        dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_index_(first_index),
        last_index_(last_index)
    { }

    basis_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int last_index
    ) : basis_(basis),
        dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_index_(0),
        last_index_(last_index)
    { }

    basis_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int first_index,
        int last_index
    ) : basis_(basis),
        dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_index_(first_index),
        last_index_(last_index)
    { }

    basis_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int first_index,
        int last_index
    ) : basis_(basis),
        dfbasis_(dfbasis),
        first_index_(first_index),
        last_index_(last_index)
    { }

    value_type begin() const;
    value_type end() const;

};

} // end namespace detail

class shell_iterator : public detail::basis_iterator<ShellData> {

  public:

    shell_iterator(const ShellBlockData& block);

    using detail::basis_iterator<ShellData>::basis_iterator;

    ShellData end() const;

};

inline shell_iterator
iter_shells_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_iterator(
      basis, dfbasis,
      shoff, shoff + basis->nshell_on_center(center) - 1
  );
}

class function_iterator : public detail::basis_iterator<BasisFunctionData> {

    int block_offset = NOT_ASSIGNED;

  public:

    using detail::basis_iterator<BasisFunctionData>::basis_iterator;

    function_iterator(const ShellData&);

    function_iterator(const ShellBlockData& block);

    BasisFunctionData begin() const;

    BasisFunctionData end() const;

};

inline function_iterator
iter_functions_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  const int bfoff = basis->shell_to_function(shoff);
  return function_iterator(
      basis, dfbasis,
      bfoff, bfoff + basis->nbasis_on_center(center) - 1
  );
}

class shell_block_iterator : public detail::basis_iterator<ShellBlockData> {

    int target_size = DEFAULT_TARGET_BLOCK_SIZE;

    // Composed using bitwise or of BlockCompositionRequirement enums
    int reqs = SameCenter;

  public:

    using detail::basis_iterator<ShellBlockData>::basis_iterator;

    const shell_block_iterator& requiring(int in_reqs) {
      reqs = in_reqs;
      return *this;
    }

    const shell_block_iterator& with_target_size(int new_target_size) {
      target_size = new_target_size;
      return *this;
    }

    ShellBlockData begin() const;
    ShellBlockData end() const;
};

inline const shell_block_iterator
iter_shell_blocks_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis = 0,
    int reqs = SameCenter
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_block_iterator(
      basis, dfbasis,
      shoff, shoff + basis->nshell_on_center(center) - 1
  ).requiring(reqs|SameCenter);
}

//############################################################################//

#define ASSERT_SHELL_BOUNDS \
  out_assert(basis, !=, 0); \
  out_assert(index, <, basis->nshell()); \
  out_assert(index, >=, 0)

struct BasisElementData {
  protected:

    BasisElementData(
        int index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis
    ) : index(index),
        basis(basis),
        dfbasis(dfbasis)
    { }

    BasisElementData() : index(NotAssigned), basis(0), dfbasis(0)  { }

  public:

    bool operator!=(const BasisElementData& other) const
    {
      return index != other.index;
    }

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

    enum { NotAssigned = -1 };

    int index = NotAssigned;

};

struct ShellData : public BasisElementData {

    ShellData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : BasisElementData(idx, basis, dfbasis)
    {
      init();
    }

    ShellData()
      : ShellData(NotAssigned, 0, 0)
    { }

    const ShellData& operator++()
    {
      ++index;
      init();
      return *this;
    }

    void set_index(int idx) {
      index = idx;
      init();
    }

    const ShellData& operator*() const { return *this; }

    int bfoff = NotAssigned;
    int nbf = NotAssigned;
    int center = NotAssigned;
    int atom_bfoff = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_nsh = NotAssigned;
    int atom_nbf = NotAssigned;
    int bfoff_in_atom = NotAssigned;
    int shoff_in_atom = NotAssigned;
    int atom_last_function = NotAssigned;
    int atom_last_shell = NotAssigned;
    int last_function = NotAssigned;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_dfnsh = NotAssigned;
    int atom_df_last_function = NotAssigned;
    int atom_df_last_shell = NotAssigned;

    operator int() { ASSERT_SHELL_BOUNDS; return index; }
    operator const int() const { ASSERT_SHELL_BOUNDS; return index; }

  private:

    void init(){
      if(index == NotAssigned || index == basis->nshell()) return;

      ASSERT_SHELL_BOUNDS;

      nbf = basis->shell(index).nfunction(); 
      bfoff = basis->shell_to_function(index); 
      center = basis->shell_to_center(index); 
      atom_shoff = basis->shell_on_center(center, 0);
      atom_bfoff = basis->shell_to_function(atom_shoff); 
      atom_nsh = basis->nshell_on_center(center); 
      atom_nbf = basis->nbasis_on_center(center); 
      shoff_in_atom = index - atom_shoff; 
      bfoff_in_atom = bfoff - atom_bfoff; 
      atom_last_function = atom_bfoff + atom_nbf - 1; 
      last_function = bfoff + nbf - 1; 
      atom_last_shell = atom_shoff + atom_nsh - 1; 
      if(dfbasis != 0) {
        atom_dfshoff = dfbasis->shell_on_center(center, 0); 
        atom_dfbfoff = dfbasis->shell_to_function(atom_dfshoff); 
        atom_dfnsh = dfbasis->nshell_on_center(center); 
        atom_dfnbf = dfbasis->nbasis_on_center(center); 
        atom_df_last_function = atom_dfbfoff + atom_dfnbf - 1; 
        atom_df_last_shell = atom_dfshoff + atom_dfnsh - 1; 
      }
    }
};

struct BasisFunctionData : public BasisElementData {

    BasisFunctionData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int block_offset = NotAssigned
    ) : BasisElementData(idx, basis, dfbasis),
        block_offset(block_offset)
    {
      init();
    }

    BasisFunctionData()
      : BasisFunctionData(NotAssigned, 0, 0)
    { }

    const BasisFunctionData& operator++()
    {
      ++index;
      init();
      return *this;
    }

    const BasisFunctionData& operator*() const { return *this; }

    int shell_index = NotAssigned;
    int shell_bfoff = NotAssigned;
    int center = NotAssigned;
    int bfoff_in_shell = NotAssigned;
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int bfoff_in_atom = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_bfoff = NotAssigned;
    int block_offset = NotAssigned;
    int bfoff_in_block = NotAssigned;

    operator int() { return index; }

  private:

    void init()
    {
      if(index == NotAssigned || index == basis->nbasis()) return;

      shell_index = basis->function_to_shell(index);
      shell_bfoff = basis->shell_to_function(shell_index);
      center = basis->shell_to_center(shell_index);
      bfoff_in_shell = index - shell_bfoff;
      atom_shoff = basis->shell_on_center(center, 0);
      atom_bfoff =  basis->shell_to_function(atom_shoff);
      bfoff_in_atom = index - atom_bfoff;
      if(dfbasis != 0){
        atom_dfshoff = dfbasis->shell_on_center(center, 0);
        atom_dfbfoff = dfbasis->shell_to_function(atom_dfshoff);
        atom_dfnbf = dfbasis->nbasis_on_center(center);
      }
      if(block_offset != NotAssigned) {
        bfoff_in_block = index - block_offset;
      }
    }
};

struct ShellBlockData {

    ShellBlockData(
        int first_shell,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int reqs = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : first_index(first_shell),
        basis(basis),
        dfbasis(dfbasis),
        target_size(target_size),
        reqs(reqs)
    {
      init();
    }

    ShellBlockData(
        const ShellBlockData& other
    ) : first_index(other.first_index),
        basis(other.basis),
        dfbasis(other.dfbasis),
        target_size(other.target_size),
        reqs(other.reqs)
    {
      init();
    }

    bool operator!=(const ShellBlockData& other) const
    {
      // Don't do last shell, since it is undefined in the case of the end() iterator
      return first_index != other.first_index or last_index != other.last_index;
    }

    const ShellBlockData& operator++()
    {
      first_index = last_index + 1;
      init();
      return *this;
    }

    void set_index(int idx) { first_index = idx; init(); }

    const ShellBlockData& operator*() const { return *this; }

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    ShellData first_shell;
    ShellData last_shell;

    int nbf;
    int bfoff;
    int last_function;
    int center = -1;
    int atom_bfoff = -1;
    int atom_shoff = -1;
    int atom_nsh = -1;
    int atom_nbf = -1;
    int bfoff_in_atom = -1;
    int shoff_in_atom = -1;
    int atom_last_function = -1;
    int atom_last_shell = -1;
    int atom_dfshoff = -1;
    int atom_dfbfoff = -1;
    int atom_dfnbf = -1;
    int atom_dfnsh = -1;
    int atom_df_last_function = -1;
    int atom_df_last_shell = -1;

  private:

    int first_index;
    int last_index;
    //std::vector<int>* shell_list = 0;
    int target_size;
    int reqs;

    void init();

};

//############################################################################//



} // end namespace sc

#include <chemistry/qc/scf/cadf_iters_impl.h>

#endif /* _chemistry_qc_scf_cadf_iters_h */
