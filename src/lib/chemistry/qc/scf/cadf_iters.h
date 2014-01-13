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

  public:

    using detail::basis_iterator<BasisFunctionData>::basis_iterator;

    function_iterator(const ShellData&);

    function_iterator(const ShellBlockData& block);

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

    BlockCompositionRequirement reqs = SameCenter;
    int target_size = DEFAULT_TARGET_BLOCK_SIZE;

  public:

    using detail::basis_iterator<ShellBlockData>::basis_iterator;

    const shell_block_iterator& requiring(BlockCompositionRequirement in_reqs) {
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


#define ShellData_USE_PROPERTIES 0
struct ShellData : public BasisElementData {

    ShellData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : BasisElementData(idx, basis, dfbasis)
        #if ShellData_USE_PROPERTIES
        ,
        bfoff(this, &ShellData::get_bfoff),
        nbf(this, &ShellData::get_nbf),
        center(this, &ShellData::get_center),
        atom_bfoff(this, &ShellData::get_atom_bfoff),
        atom_shoff(this, &ShellData::get_atom_shoff),
        atom_nsh(this, &ShellData::get_atom_nsh),
        atom_nbf(this, &ShellData::get_atom_nbf),
        bfoff_in_atom(this, &ShellData::get_bfoff_in_atom),
        shoff_in_atom(this, &ShellData::get_shoff_in_atom),
        atom_last_function(this, &ShellData::get_atom_last_function),
        atom_last_shell(this, &ShellData::get_atom_last_shell),
        last_function(this, &ShellData::get_last_function),
        atom_dfshoff(this, &ShellData::get_atom_dfshoff),
        atom_dfbfoff(this, &ShellData::get_atom_dfbfoff),
        atom_dfnbf(this, &ShellData::get_atom_dfnbf),
        atom_dfnsh(this, &ShellData::get_atom_dfnsh),
        atom_df_last_function(this, &ShellData::get_atom_df_last_function),
        atom_df_last_shell(this, &ShellData::get_atom_df_last_shell)
        #endif
    {
      #if !ShellData_USE_PROPERTIES
      init();
      #endif
    }

    ShellData()
      : ShellData(NotAssigned, 0, 0)
    { }

    ShellData(const ShellData& other)
      : ShellData(other.index, other.basis, other.dfbasis)
    { }

    bool operator!=(const ShellData& other) const
    {
      return index != other.index;
    }

    const ShellData& operator++()
    {
      ++index;
      #if !ShellData_USE_PROPERTIES
      init();
      #endif
      return *this;
    }

    void set_index(int idx) {
      index = idx;
      #if !ShellData_USE_PROPERTIES
      init();
      #endif
    }

    const ShellData& operator*() const { return *this; }

    ShellData&
    operator=(const ShellData& other)
    {
      index = other.index;
      basis = other.basis;
      dfbasis = other.dfbasis;
      #if ShellData_USE_PROPERTIES
      bfoff.set_target(this);
      nbf.set_target(this);
      center.set_target(this);
      atom_bfoff.set_target(this);
      atom_shoff.set_target(this);
      atom_nsh.set_target(this);
      atom_nbf.set_target(this);
      bfoff_in_atom.set_target(this);
      shoff_in_atom.set_target(this);
      atom_last_function.set_target(this);
      atom_last_shell.set_target(this);
      last_function.set_target(this);
      atom_dfshoff.set_target(this);
      atom_dfbfoff.set_target(this);
      atom_dfnbf.set_target(this);
      atom_dfnsh.set_target(this);
      atom_df_last_function.set_target(this);
      atom_df_last_shell.set_target(this);
      #else
      init();
      #endif
      return *this;
    }

    #if ShellData_USE_PROPERTIES
    template<typename t1, typename t2> using _property = property<t1, t2>;
    #else
    template<typename t1, typename t2> using _property = t2;
    #endif
    _property<const ShellData, int> bfoff;
    _property<const ShellData, int> nbf;
    _property<const ShellData, int> center;
    _property<const ShellData, int> atom_bfoff;
    _property<const ShellData, int> atom_shoff;
    _property<const ShellData, int> atom_nsh;
    _property<const ShellData, int> atom_nbf;
    _property<const ShellData, int> bfoff_in_atom;
    _property<const ShellData, int> shoff_in_atom;
    _property<const ShellData, int> atom_last_function;
    _property<const ShellData, int> atom_last_shell;
    _property<const ShellData, int> last_function;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    _property<const ShellData, int> atom_dfshoff;
    _property<const ShellData, int> atom_dfbfoff;
    _property<const ShellData, int> atom_dfnbf;
    _property<const ShellData, int> atom_dfnsh;
    _property<const ShellData, int> atom_df_last_function;
    _property<const ShellData, int> atom_df_last_shell;

    operator int() { ASSERT_SHELL_BOUNDS; return index; }
    operator const int() const { ASSERT_SHELL_BOUNDS; return index; }

  private:

    #if ShellData_USE_PROPERTIES
    int get_nbf() const { ASSERT_SHELL_BOUNDS; return basis->shell(index).nfunction(); }
    int get_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(index); }
    int get_center() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_center(index); }
    int get_atom_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(atom_shoff); }
    int get_atom_shoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_on_center(center, 0); }
    int get_atom_nsh() const { ASSERT_SHELL_BOUNDS; return basis->nshell_on_center(center); }
    int get_atom_nbf() const { ASSERT_SHELL_BOUNDS; return basis->nbasis_on_center(center); }
    int get_shoff_in_atom() const { ASSERT_SHELL_BOUNDS; return index - atom_shoff; }
    int get_bfoff_in_atom() const { ASSERT_SHELL_BOUNDS; return bfoff - atom_bfoff; }
    int get_atom_last_function() const { ASSERT_SHELL_BOUNDS; return atom_bfoff + atom_nbf - 1; }
    int get_last_function() const { ASSERT_SHELL_BOUNDS; return bfoff + nbf - 1; }
    int get_atom_last_shell() const { ASSERT_SHELL_BOUNDS; return atom_shoff + atom_nsh - 1; }
    int get_atom_dfshoff() const { assert(dfbasis != 0); return dfbasis->shell_on_center(center, 0); }
    int get_atom_dfbfoff() const { assert(dfbasis != 0); return dfbasis->shell_to_function(atom_dfshoff); }
    int get_atom_dfnsh() const { assert(dfbasis != 0); return dfbasis->nshell_on_center(center); }
    int get_atom_dfnbf() const { assert(dfbasis != 0); return dfbasis->nbasis_on_center(center); }
    int get_atom_df_last_function() const { assert(dfbasis != 0); return atom_dfbfoff + atom_dfnbf - 1; }
    int get_atom_df_last_shell() const { assert(dfbasis != 0); return atom_dfshoff + atom_dfnsh - 1; }
    int assert_not_initialized() const { assert(false && "ShellData object not initialized"); return -1; }
    #else
    void init(){
      if(index == NotAssigned || index == basis->nshell()) return;

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
    #endif
};


#define BFD_USE_PROPERTIES 0
struct BasisFunctionData {

    BasisFunctionData()
      : index(-1),
        basis(0),
        dfbasis(0)
        #if BFD_USE_PROPERTIES
        ,
        shell_index(this, &BasisFunctionData::assert_not_initialized),
        shell_bfoff(this, &BasisFunctionData::assert_not_initialized),
        center(this, &BasisFunctionData::assert_not_initialized),
        bfoff_in_shell(this, &BasisFunctionData::assert_not_initialized),
        atom_dfshoff(this, &BasisFunctionData::assert_not_initialized),
        atom_dfbfoff(this, &BasisFunctionData::assert_not_initialized),
        atom_dfnbf(this, &BasisFunctionData::assert_not_initialized),
        bfoff_in_atom(this, &BasisFunctionData::assert_not_initialized),
        atom_bfoff(this, &BasisFunctionData::assert_not_initialized),
        atom_shoff(this, &BasisFunctionData::assert_not_initialized)
        //shell_nbf(this, &BasisFunctionData::assert_not_initialized),
        //atom_nsh(this, &BasisFunctionData::assert_not_initialized),
        //atom_nbf(this, &BasisFunctionData::assert_not_initialized),
        //shoff_in_atom(this, &BasisFunctionData::assert_not_initialized),
        //atom_last_function(this, &BasisFunctionData::assert_not_initialized),
        //atom_last_shell(this, &BasisFunctionData::assert_not_initialized),
        //atom_dfnsh(this, &BasisFunctionData::assert_not_initialized),
        //atom_df_last_function(this, &BasisFunctionData::assert_not_initialized),
        //atom_df_last_shell(this, &BasisFunctionData::assert_not_initialized)
        #endif
    {

    }

    BasisFunctionData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis
    ) : index(idx),
        basis(basis),
        dfbasis(dfbasis)
        #if BFD_USE_PROPERTIES
        ,
        shell_index(this, &BasisFunctionData::get_shell_index),
        shell_bfoff(this, &BasisFunctionData::get_shell_bfoff),
        center(this, &BasisFunctionData::get_center),
        bfoff_in_shell(this, &BasisFunctionData::get_bfoff_in_shell),
        atom_dfshoff(this, &BasisFunctionData::get_atom_dfshoff),
        atom_dfbfoff(this, &BasisFunctionData::get_atom_dfbfoff),
        atom_dfnbf(this, &BasisFunctionData::get_atom_dfnbf),
        bfoff_in_atom(this, &BasisFunctionData::get_bfoff_in_atom),
        atom_bfoff(this, &BasisFunctionData::get_atom_bfoff),
        atom_shoff(this, &BasisFunctionData::get_atom_shoff)
        //shell_nbf(this, &BasisFunctionData::get_shell_nbf),
        //atom_nsh(this, &BasisFunctionData::get_atom_nsh),
        //atom_nbf(this, &BasisFunctionData::get_atom_nbf),
        //shoff_in_atom(this, &BasisFunctionData::get_shoff_in_atom),
        //atom_last_function(this, &BasisFunctionData::get_atom_last_function),
        //atom_last_shell(this, &BasisFunctionData::get_atom_last_shell),
        //atom_dfnsh(this, &BasisFunctionData::get_atom_dfnsh),
        //atom_df_last_function(this, &BasisFunctionData::get_atom_df_last_function),
        //atom_df_last_shell(this, &BasisFunctionData::get_atom_df_last_shell)
        #endif
    {
      #if BFD_USE_PROPERTIES == 0
      init();
      #endif
    }

    BasisFunctionData(
        const BasisFunctionData& other
    ) : index(other.index),
        basis(other.basis),
        dfbasis(other.dfbasis)
        #if BFD_USE_PROPERTIES
        ,
        shell_index(this, &BasisFunctionData::get_shell_index),
        shell_bfoff(this, &BasisFunctionData::get_shell_bfoff),
        center(this, &BasisFunctionData::get_center),
        bfoff_in_shell(this, &BasisFunctionData::get_bfoff_in_shell),
        atom_dfshoff(this, &BasisFunctionData::get_atom_dfshoff),
        atom_dfbfoff(this, &BasisFunctionData::get_atom_dfbfoff),
        atom_dfnbf(this, &BasisFunctionData::get_atom_dfnbf),
        bfoff_in_atom(this, &BasisFunctionData::get_bfoff_in_atom),
        atom_bfoff(this, &BasisFunctionData::get_atom_bfoff),
        atom_shoff(this, &BasisFunctionData::get_atom_shoff)
        //shell_nbf(this, &BasisFunctionData::get_shell_nbf),
        //atom_nsh(this, &BasisFunctionData::get_atom_nsh),
        //atom_nbf(this, &BasisFunctionData::get_atom_nbf),
        //shoff_in_atom(this, &BasisFunctionData::get_shoff_in_atom),
        //atom_last_function(this, &BasisFunctionData::get_atom_last_function),
        //atom_last_shell(this, &BasisFunctionData::get_atom_last_shell),
        //atom_dfnsh(this, &BasisFunctionData::get_atom_dfnsh),
        //atom_df_last_function(this, &BasisFunctionData::get_atom_df_last_function),
        //atom_df_last_shell(this, &BasisFunctionData::get_atom_df_last_shell)
        #endif
    {
      #if BFD_USE_PROPERTIES == 0
      init();
      #endif
    }

    bool operator!=(const BasisFunctionData& other) const
    {
      return index != other.index;
    }

    const BasisFunctionData& operator++()
    {
      ++index;
      #if BFD_USE_PROPERTIES == 0
      init();
      #endif
      return *this;
    }

    const BasisFunctionData& operator*() const { return *this; }

    int index;

    #if BFD_USE_PROPERTIES
    property<BasisFunctionData, int> shell_index;
    property<BasisFunctionData, int> shell_bfoff;
    property<BasisFunctionData, int> center;
    property<BasisFunctionData, int> bfoff_in_shell;
    property<BasisFunctionData, int> atom_dfshoff;
    property<BasisFunctionData, int> atom_dfbfoff;
    property<BasisFunctionData, int> atom_dfnbf;
    property<BasisFunctionData, int> bfoff_in_atom;
    property<BasisFunctionData, int> atom_bfoff;
    property<BasisFunctionData, int> atom_shoff;
    //property<BasisFunctionData, int> shell_nbf;
    //property<BasisFunctionData, int> atom_nsh;
    //property<BasisFunctionData, int> atom_nbf;
    //property<BasisFunctionData, int> shoff_in_atom;
    //property<BasisFunctionData, int> atom_last_function;
    //property<BasisFunctionData, int> atom_last_shell;
    //property<BasisFunctionData, int> atom_dfnsh;
    //property<BasisFunctionData, int> atom_df_last_function;
    //property<BasisFunctionData, int> atom_df_last_shell;
    #else
    int shell_index = -1;
    int shell_bfoff = -1;
    int center = -1;
    int bfoff_in_shell = -1;
    int atom_dfshoff = -1;
    int atom_dfbfoff = -1;
    int atom_dfnbf = -1;
    int bfoff_in_atom = -1;
    int atom_shoff = -1;
    int atom_bfoff = -1;
    #endif

    operator int() { return index; }

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

  private:

    #if not BFD_USE_PROPERTIES
    void init()
    {
      //std::cout << "index = " << index << ", " << "nbasis() = " << basis->nbasis() << std::endl;
      if(index < basis->nbasis() && index >= 0){
        shell_index = get_shell_index();
        //std::cout << "  shell_index = " << shell_index << std::endl;
        shell_bfoff = get_shell_bfoff();
        center = get_center();
        bfoff_in_shell = get_bfoff_in_shell();
        atom_shoff = get_atom_shoff();
        atom_bfoff = get_atom_bfoff();
        bfoff_in_atom = get_bfoff_in_atom();
        if(dfbasis != 0){
          atom_dfshoff = get_atom_dfshoff();
          atom_dfbfoff = get_atom_dfbfoff();
          atom_dfnbf = get_atom_dfnbf();
        }
      }
    }
    #endif

    int get_shell_index() const { return basis->function_to_shell(index); }
    int get_shell_bfoff() const { return basis->shell_to_function(shell_index); }
    int get_center() const { return basis->shell_to_center(shell_index); }
    int get_bfoff_in_shell() const { return index - shell_bfoff; }
    int get_atom_dfshoff() const { assert(dfbasis != 0); return dfbasis->shell_on_center(center, 0); }
    int get_atom_dfbfoff() const { assert(dfbasis != 0); return dfbasis->shell_to_function(atom_dfshoff); }
    int get_atom_dfnbf() const { assert(dfbasis != 0); return dfbasis->nbasis_on_center(center); }
    int assert_not_initialized() const { assert(false && "ShellData object not initialized"); return -1; }
    int get_bfoff_in_atom() const { return index - atom_bfoff; }
    int get_atom_bfoff() const { return basis->shell_to_function(atom_shoff); }
    int get_atom_shoff() const { return basis->shell_on_center(center, 0); }
    //int get_shell_nbf() const { return basis->shell(shell_index).nfunction(); }
    //int get_atom_nsh() const { return basis->nshell_on_center(center); }
    //int get_atom_nbf() const { return basis->nbasis_on_center(center); }
    //int get_shoff_in_atom() const { return shell_index - atom_shoff; }
    //int get_atom_last_function() const { return atom_bfoff + atom_nbf - 1; }
    //int get_atom_last_shell() const { return atom_shoff + atom_nsh - 1; }
    //int get_atom_dfnsh() const { assert(dfbasis.nonnull()); return dfbasis->nshell_on_center(center); }
    //int get_atom_df_last_function() const { assert(dfbasis.nonnull()); return atom_dfbfoff + atom_dfnbf - 1; }
    //int get_atom_df_last_shell() const { assert(dfbasis.nonnull()); return atom_dfshoff + atom_dfnsh - 1; }
};

struct ShellBlockData {


    ShellBlockData(
        int first_shell,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        BlockCompositionRequirement reqs = SameCenter,
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
    BlockCompositionRequirement reqs;

    void init();

};

//############################################################################//



} // end namespace sc

#include <chemistry/qc/scf/cadf_iters_impl.h>

#endif /* _chemistry_qc_scf_cadf_iters_h */
