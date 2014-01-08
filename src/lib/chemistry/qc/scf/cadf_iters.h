//
// cadf_iters.h
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Dec 18, 2013
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

#define out_assert(a, op, b) assert(a op b || ((std::cout << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))

namespace sc {

// Forward declarations
class ShellData;
class BasisFunctionData;

class shell_iterator {
    GaussianBasisSet* basis_;
    GaussianBasisSet* dfbasis_;
    int first_shell_;
    int last_shell_;
  public:
    shell_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0), first_shell_(0), last_shell_(-1)
    { }

    shell_iterator(
        const Ref<GaussianBasisSet>& basis,
        int last_shell,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0), first_shell_(0), last_shell_(last_shell)
    { }

    shell_iterator(
        const Ref<GaussianBasisSet>& basis,
        int first_shell,
        int last_shell,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0), first_shell_(first_shell), last_shell_(last_shell)
    { }

    shell_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int last_shell
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0), first_shell_(0), last_shell_(last_shell)
    { }

    shell_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int first_shell,
        int last_shell
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0), first_shell_(first_shell), last_shell_(last_shell)
    { }

    ShellData begin() const;
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

class function_iterator {

    GaussianBasisSet* basis_;
    GaussianBasisSet* dfbasis_;
    int first_function_;
    int last_function_;

  public:
    function_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_function_(0), last_function_(-1)
    { }

    function_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis),
        first_function_(0), last_function_(-1)
    { }

    function_iterator(
        const Ref<GaussianBasisSet>& basis,
        int last_function,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_function_(0), last_function_(last_function)
    { }

    function_iterator(
        const Ref<GaussianBasisSet>& basis,
        int first_function,
        int last_function,
        const Ref<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_function_(first_function), last_function_(last_function)
    { }

    function_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int last_function
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_function_(0), last_function_(last_function)
    { }

    function_iterator(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int first_function,
        int last_function
    ) : basis_(basis), dfbasis_(dfbasis.nonnull() ? dfbasis.pointer() : 0),
        first_function_(first_function), last_function_(last_function)
    { }

    function_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int first_function,
        int last_function
    ) : basis_(basis), dfbasis_(dfbasis),
        first_function_(first_function), last_function_(last_function)
    { }

    function_iterator(const ShellData&);

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


//============================================================================//

#define ASSERT_SHELL_BOUNDS \
  out_assert(basis, !=, 0); \
  out_assert(index, <, basis->nshell()); \
  out_assert(index, >=, 0)

struct ShellData {

    ShellData()
      : index(-1),
        basis(0),
        dfbasis(0),
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
    {

    }

    ShellData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : index(idx),
        basis(basis),
        dfbasis(dfbasis),
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
    {

    }

    ShellData(const ShellData& other)
      : index(other.index),
        basis(other.basis),
        dfbasis(other.dfbasis),
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
    {

    }

    bool operator!=(const ShellData& other) const
    {
      return index != other.index;
    }

    const ShellData& operator++()
    {
      ++index;
      return *this;
    }

    void set_index(int idx) { index = idx; }

    const ShellData& operator*() const { return *this; }

    ShellData&
    operator=(const ShellData& other)
    {
      index = other.index;
      basis = other.basis;
      dfbasis = other.dfbasis;
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
      return *this;
    }

    int index;
    property<const ShellData, int> bfoff;
    property<const ShellData, int> nbf;
    property<const ShellData, int> center;
    property<const ShellData, int> atom_bfoff;
    property<const ShellData, int> atom_shoff;
    property<const ShellData, int> atom_nsh;
    property<const ShellData, int> atom_nbf;
    property<const ShellData, int> bfoff_in_atom;
    property<const ShellData, int> shoff_in_atom;
    property<const ShellData, int> atom_last_function;
    property<const ShellData, int> atom_last_shell;
    property<const ShellData, int> last_function;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    property<const ShellData, int> atom_dfshoff;
    property<const ShellData, int> atom_dfbfoff;
    property<const ShellData, int> atom_dfnbf;
    property<const ShellData, int> atom_dfnsh;
    property<const ShellData, int> atom_df_last_function;
    property<const ShellData, int> atom_df_last_shell;

    operator int() { ASSERT_SHELL_BOUNDS; return index; }

    //const Ref<GaussianBasisSet>& basis;
    //const Ref<GaussianBasisSet>& dfbasis;
    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

  private:
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

//============================================================================//

inline
function_iterator::function_iterator(const ShellData& ish)
  : basis_(ish.basis), dfbasis_(ish.dfbasis), first_function_(ish.bfoff), last_function_(ish.last_function)
{
  /*
  assert(ish.basis != 0);
  std::cout << "ish.basis->nbasis() = " << ish.basis->nbasis() << std::endl;
  */
}

inline BasisFunctionData
function_iterator::begin() const
{
  return BasisFunctionData(first_function_, basis_, dfbasis_);
}

inline BasisFunctionData
function_iterator::end() const
{
  return BasisFunctionData(last_function_ == -1 ? basis_->nbasis() : last_function_ + 1, basis_, dfbasis_);
}

//============================================================================//
// ShellData

inline ShellData
shell_data(
    Ref<GaussianBasisSet> basis,
    int ish, Ref<GaussianBasisSet> dfbasis
)
{
  return ShellData(ish, basis, dfbasis);
}

/*inline ShellData&
ShellData::operator=(const ShellData& other)
{
  index = other.index;
  basis = other.basis;
  dfbasis = other.dfbasis;
  bfoff.set_target_getter_setter(this, &ShellData::get_bfoff);
  nbf.set_target_getter_setter(this, &ShellData::get_nbf);
  center.set_target_getter_setter(this, &ShellData::get_center);
  atom_bfoff.set_target_getter_setter(this, &ShellData::get_atom_bfoff);
  atom_shoff.set_target_getter_setter(this, &ShellData::get_atom_shoff);
  atom_nsh.set_target_getter_setter(this, &ShellData::get_atom_nsh);
  atom_nbf.set_target_getter_setter(this, &ShellData::get_atom_nbf);
  bfoff_in_atom.set_target_getter_setter(this, &ShellData::get_bfoff_in_atom);
  shoff_in_atom.set_target_getter_setter(this, &ShellData::get_shoff_in_atom);
  atom_last_function.set_target_getter_setter(this, &ShellData::get_atom_last_function);
  atom_last_shell.set_target_getter_setter(this, &ShellData::get_atom_last_shell);
  last_function.set_target_getter_setter(this, &ShellData::get_last_function);
  atom_dfshoff.set_target_getter_setter(this, &ShellData::get_atom_dfshoff);
  atom_dfbfoff.set_target_getter_setter(this, &ShellData::get_atom_dfbfoff);
  atom_dfnbf.set_target_getter_setter(this, &ShellData::get_atom_dfnbf);
  atom_dfnsh.set_target_getter_setter(this, &ShellData::get_atom_dfnsh);
  atom_df_last_function.set_target_getter_setter(this, &ShellData::get_atom_df_last_function);
  atom_df_last_shell.set_target_getter_setter(this, &ShellData::get_atom_df_last_shell);
  return *this;
}
*/

inline const BasisFunctionData
function_data(
    const Ref<GaussianBasisSet>& basis,
    int ish, const Ref<GaussianBasisSet>& dfbasis
)
{
  return BasisFunctionData(ish, basis, dfbasis);
}

inline ShellData
shell_iterator::begin() const
{
  return ShellData(first_shell_, basis_, dfbasis_);
}

inline ShellData
shell_iterator::end() const
{
  return ShellData(last_shell_ == -1 ? basis_->nshell() : last_shell_ + 1, basis_, dfbasis_);
}

template <typename Iterator>
struct shell_iter_arbitrary {

    shell_iter_arbitrary(
        Iterator iter,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : iter_(iter),
        sd((const int)*iter, basis, dfbasis)
    { }

    const shell_iter_arbitrary& operator++() {
      ++iter_;
      sd.set_index((int)(*iter_));
      return *this;
    }

    bool operator!=(const shell_iter_arbitrary& other) const {
      return iter_ != other.iter_;
    }

    const ShellData& operator*() const { return sd; }

  private:
    ShellData sd;
    Iterator iter_;

    BOOST_STATIC_ASSERT(
        (std::is_convertible<decltype(iter_.operator*()), int>::value)
    );
};

template <typename Iterable>
struct shell_iter_arbitrary_wrapper {
    shell_iter_arbitrary_wrapper(
        const Iterable& iterable,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis
    ) : iterable(iterable), basis(basis), dfbasis(dfbasis)
    {

    }

    shell_iter_arbitrary<typename Iterable::const_iterator>
    begin() const {
      return shell_iter_arbitrary<typename Iterable::const_iterator>(
          iterable.begin(), basis, dfbasis);
    }

    shell_iter_arbitrary<typename Iterable::const_iterator>
    end() const {
      return shell_iter_arbitrary<typename Iterable::const_iterator>(
          iterable.end(), basis, dfbasis);
    }

  private:
    const Iterable& iterable;
    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
};



//============================================================================//

} // end namespace sc



#endif /* _chemistry_qc_scf_cadf_iters_h */
