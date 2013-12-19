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

namespace sc {

// Forward declarations
class ShellData;
class BasisFunctionData;

class shell_iter_wrapper {
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;
    int first_shell_;
    int last_shell_;
  public:
    shell_iter_wrapper(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0,
        int first_shell = 0,
        int last_shell = -1
    );
    ShellIter begin() const;
    ShellIter end() const;
};

class function_iter_wrapper {
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;
    int first_function_;
    int last_function_;
  public:
    function_iter_wrapper(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0,
        int first_function = 0,
        int last_function = -1
    );
    BasisFunctionIter begin() const;
    BasisFunctionIter end() const;
};

class ShellIter {

  public:

    ShellIter() = delete;

    ShellIter(const Ref<GaussianBasisSet>& basis, int position);
    ShellIter(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int position);

    bool operator!=(const ShellIter& other) const;

    const ShellIter& operator++();

    const ShellData operator*() const;

  private:

    int pos_;
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;

};

class BasisFunctionIter {

  public:

    BasisFunctionIter() = delete;

    BasisFunctionIter(const Ref<GaussianBasisSet>& basis, int position);
    BasisFunctionIter(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int position);

    bool operator!=(const BasisFunctionIter& other) const;

    const BasisFunctionIter& operator++();

    const BasisFunctionData operator*() const;

  private:

    int pos_;
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;

};

//============================================================================//
// function_iterator_wrapper

inline
function_iter_wrapper::function_iter_wrapper(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_function,
    int last_function
) : basis_(basis),
    dfbasis_(dfbasis),
    first_function_(first_function),
    last_function_(last_function)
{

}

inline BasisFunctionIter
function_iter_wrapper::begin() const
{
  if(dfbasis_.nonnull()){
    BasisFunctionIter rv(basis_, dfbasis_, first_function_);
    return rv;
  }
  else{
    BasisFunctionIter rv(basis_, first_function_);
    return rv;
  }
}

inline BasisFunctionIter
function_iter_wrapper::end() const
{
  if(dfbasis_.nonnull()){
    BasisFunctionIter rv(basis_, dfbasis_, last_function_ == -1 ? basis_->nbasis() : last_function_ + 1);
    return rv;
  }
  else{
    BasisFunctionIter rv(basis_, last_function_ == -1 ? basis_->nbasis() : last_function_ + 1);
    return rv;
  }
}

//============================================================================//
// shell_iterator

inline
const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis = 0)
{
  return shell_iter_wrapper(basis, dfbasis, 0, -1);
}

inline
const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, int last_shell)
{
  return shell_iter_wrapper(basis, 0, 0, last_shell);
}

inline
const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, int first_shell, int last_shell)
{
  return shell_iter_wrapper(basis, 0, first_shell, last_shell);
}

inline
const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int last_shell)
{
  return shell_iter_wrapper(basis, dfbasis, 0, last_shell);
}

inline
const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int first_shell, int last_shell)
{
  return shell_iter_wrapper(basis, dfbasis, first_shell, last_shell);
}

//============================================================================//
// function_iterator

inline
const function_iter_wrapper
function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis = 0
)
{
  return function_iter_wrapper(basis, dfbasis, 0, -1);
}

inline
const function_iter_wrapper
function_iterator(
    const Ref<GaussianBasisSet>& basis,
    int last_function,
    const Ref<GaussianBasisSet>& dfbasis = 0
)
{
  return function_iter_wrapper(basis, dfbasis, 0, last_function);
}

inline
const function_iter_wrapper
function_iterator(
    const Ref<GaussianBasisSet>& basis,
    int first_function,
    int last_function,
    const Ref<GaussianBasisSet>& dfbasis = 0
)
{
  return function_iter_wrapper(basis, dfbasis, first_function, last_function);
}

inline
const function_iter_wrapper
function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_function,
    int last_function
)
{
  return function_iter_wrapper(basis, dfbasis, first_function, last_function);
}

inline
const function_iter_wrapper
function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int last_function
)
{
  return function_iter_wrapper(basis, dfbasis, 0, last_function);
}

#define ASSERT_SHELL_BOUNDS \
  assert(basis.nonnull() && "null basis"); \
  assert(index < basis->nshell() && "index larger than nshell()"); \
  assert(index >= 0 && "index less than 0")

struct ShellData {

    ShellData();

    ShellData(
        int idx,
        Ref<GaussianBasisSet> basis,
        Ref<GaussianBasisSet> dfbasis = 0
    );

    ShellData(const ShellData&);

    ShellData& operator=(const ShellData& other);

    int index;
    property<ShellData, int> bfoff;
    property<ShellData, int> nbf;
    property<ShellData, int> center;
    property<ShellData, int> atom_bfoff;
    property<ShellData, int> atom_shoff;
    property<ShellData, int> atom_nsh;
    property<ShellData, int> atom_nbf;
    property<ShellData, int> bfoff_in_atom;
    property<ShellData, int> shoff_in_atom;
    property<ShellData, int> atom_last_function;
    property<ShellData, int> atom_last_shell;
    property<ShellData, int> last_function;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    property<ShellData, int> atom_dfshoff;
    property<ShellData, int> atom_dfbfoff;
    property<ShellData, int> atom_dfnbf;
    property<ShellData, int> atom_dfnsh;
    property<ShellData, int> atom_df_last_function;
    property<ShellData, int> atom_df_last_shell;

    operator int() { ASSERT_SHELL_BOUNDS; return index; }

    Ref<GaussianBasisSet> basis;
    Ref<GaussianBasisSet> dfbasis;

  private:
    inline int get_nbf() const { ASSERT_SHELL_BOUNDS; return basis->shell(index).nfunction(); }
    inline int get_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(index); }
    inline int get_center() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_center(index); }
    inline int get_atom_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(atom_shoff); }
    inline int get_atom_shoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_on_center(center, 0); }
    inline int get_atom_nsh() const { ASSERT_SHELL_BOUNDS; return basis->nshell_on_center(center); }
    inline int get_atom_nbf() const { ASSERT_SHELL_BOUNDS; return basis->nbasis_on_center(center); }
    inline int get_shoff_in_atom() const { ASSERT_SHELL_BOUNDS; return index - atom_shoff; }
    inline int get_bfoff_in_atom() const { ASSERT_SHELL_BOUNDS; return bfoff - atom_bfoff; }
    inline int get_atom_last_function() const { ASSERT_SHELL_BOUNDS; return atom_bfoff + atom_nbf - 1; }
    inline int get_last_function() const { ASSERT_SHELL_BOUNDS; return bfoff + nbf - 1; }
    inline int get_atom_last_shell() const { ASSERT_SHELL_BOUNDS; return atom_shoff + atom_nsh - 1; }
    inline int get_atom_dfshoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_on_center(center, 0); }
    inline int get_atom_dfbfoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_to_function(atom_dfshoff); }
    inline int get_atom_dfnsh() const { assert(dfbasis.nonnull()); return dfbasis->nshell_on_center(center); }
    inline int get_atom_dfnbf() const { assert(dfbasis.nonnull()); return dfbasis->nbasis_on_center(center); }
    inline int get_atom_df_last_function() const { assert(dfbasis.nonnull()); return atom_dfbfoff + atom_dfnbf - 1; }
    inline int get_atom_df_last_shell() const { assert(dfbasis.nonnull()); return atom_dfshoff + atom_dfnsh - 1; }
    inline int assert_not_initialized() const { assert(false && "ShellData object not initialized"); return -1; }

};

struct BasisFunctionData {

    BasisFunctionData();

    BasisFunctionData(
        int idx,
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0
    );

    BasisFunctionData(const BasisFunctionData&);

    BasisFunctionData& operator=(const BasisFunctionData& other);

    int index;
    property<BasisFunctionData, int> shell_index;
    property<BasisFunctionData, int> shell_nbf;
    property<BasisFunctionData, int> shell_bfoff;
    property<BasisFunctionData, int> center;
    property<BasisFunctionData, int> atom_bfoff;
    property<BasisFunctionData, int> atom_shoff;
    property<BasisFunctionData, int> atom_nsh;
    property<BasisFunctionData, int> atom_nbf;
    property<BasisFunctionData, int> bfoff_in_atom;
    property<BasisFunctionData, int> bfoff_in_shell;
    property<BasisFunctionData, int> shoff_in_atom;
    property<BasisFunctionData, int> atom_last_function;
    property<BasisFunctionData, int> atom_last_shell;
    property<BasisFunctionData, int> atom_dfshoff;
    property<BasisFunctionData, int> atom_dfbfoff;
    property<BasisFunctionData, int> atom_dfnbf;
    property<BasisFunctionData, int> atom_dfnsh;
    property<BasisFunctionData, int> atom_df_last_function;
    property<BasisFunctionData, int> atom_df_last_shell;

    operator int() { return index; }

    Ref<GaussianBasisSet> basis;
    Ref<GaussianBasisSet> dfbasis;

  private:
    inline int get_shell_index() const { return basis->function_to_shell(index); }
    inline int get_shell_nbf() const { return basis->shell(shell_index).nfunction(); }
    inline int get_shell_bfoff() const { return basis->shell_to_function(shell_index); }
    inline int get_center() const { return basis->shell_to_center(shell_index); }
    inline int get_atom_bfoff() const { return basis->shell_to_function(atom_shoff); }
    inline int get_atom_shoff() const { return basis->shell_on_center(center, 0); }
    inline int get_atom_nsh() const { return basis->nshell_on_center(center); }
    inline int get_atom_nbf() const { return basis->nbasis_on_center(center); }
    inline int get_bfoff_in_atom() const { return index - atom_bfoff; }
    inline int get_bfoff_in_shell() const { return index - shell_bfoff; }
    inline int get_shoff_in_atom() const { return shell_index - atom_shoff; }
    inline int get_atom_last_function() const { return atom_bfoff + atom_nbf - 1; }
    inline int get_atom_last_shell() const { return atom_shoff + atom_nsh - 1; }
    inline int get_atom_dfshoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_on_center(center, 0); }
    inline int get_atom_dfbfoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_to_function(atom_dfshoff); }
    inline int get_atom_dfnbf() const { assert(dfbasis.nonnull()); return dfbasis->nbasis_on_center(center); }
    inline int get_atom_dfnsh() const { assert(dfbasis.nonnull()); return dfbasis->nshell_on_center(center); }
    inline int get_atom_df_last_function() const { assert(dfbasis.nonnull()); return atom_dfbfoff + atom_dfnbf - 1; }
    inline int get_atom_df_last_shell() const { assert(dfbasis.nonnull()); return atom_dfshoff + atom_dfnsh - 1; }
    inline int assert_not_initialized() const { assert(false && "ShellData object not initialized"); return -1; }
};

inline const function_iter_wrapper
function_iterator(const ShellData& ish)
{
  return function_iter_wrapper(ish.basis, ish.dfbasis, ish.bfoff, ish.last_function);
}

const shell_iter_wrapper
iter_shells_on_center(const Ref<GaussianBasisSet>& basis, int center, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
iter_functions_on_center(const Ref<GaussianBasisSet>& basis, int center, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
iter_functions_in_shell(const Ref<GaussianBasisSet>& basis, int shell, const Ref<GaussianBasisSet>& dfbasis = 0);

ShellData shell_data(Ref<GaussianBasisSet> basis, int ish, Ref<GaussianBasisSet> dfbasis = 0);

const BasisFunctionData function_data(const Ref<GaussianBasisSet>& basis, int ish, const Ref<GaussianBasisSet>& dfbasis = 0);

//============================================================================//
// ShellIter

inline
ShellIter::ShellIter(
    const Ref<GaussianBasisSet>& basis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(0)
{

}

inline
ShellIter::ShellIter(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(dfbasis)
{

}

inline bool
ShellIter::operator!=(const ShellIter& other) const
{
  return pos_ != other.pos_;
  // This doesn't seem to work. Maybe there's a mistake?  It's not that important
  //    or not basis_->equiv(other.basis_)
  //    or (dfbasis_.nonnull() and not other.basis_.nonnull())
  //    or (dfbasis_.null() and not other.basis_.null())
  //    or (dfbasis_.nonnull() and not dfbasis_->equiv(other.dfbasis_));
}

inline const ShellIter&
ShellIter::operator++()
{
  ++pos_;
  return *this;
}

inline const ShellData
ShellIter::operator*() const
{
  if(dfbasis_.nonnull()) {
    ShellData rv(pos_, basis_, dfbasis_);
    return rv;
  }
  else {
    ShellData rv(pos_, basis_);
    return rv;
  }

}

//============================================================================//
// BasisFunctionIter

inline
BasisFunctionIter::BasisFunctionIter(
    const Ref<GaussianBasisSet>& basis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(0)
{

}

inline
BasisFunctionIter::BasisFunctionIter(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(dfbasis)
{

}

inline bool
BasisFunctionIter::operator!=(const BasisFunctionIter& other) const
{
  return pos_ != other.pos_;
  // This doesn't seem to work. Maybe there's a mistake?  It's not that important
  //    or not basis_->equiv(other.basis_)
  //    or (dfbasis_.nonnull() and not other.basis_.nonnull())
  //    or (dfbasis_.null() and not other.basis_.null())
  //    or (dfbasis_.nonnull() and not dfbasis_->equiv(other.dfbasis_));
}

inline const BasisFunctionIter&
BasisFunctionIter::operator++()
{
  ++pos_;
  return *this;
}

inline const BasisFunctionData
BasisFunctionIter::operator*() const
{
  if(dfbasis_.nonnull()) {
    BasisFunctionData rv(pos_, basis_, dfbasis_);
    return rv;
  }
  else {
    BasisFunctionData rv(pos_, basis_);
    return rv;
  }
}

//============================================================================//
// ShellData

inline
ShellData::ShellData()
  : index(-1),
    basis(0),
    dfbasis(0),
    bfoff(this, &ShellData::assert_not_initialized),
    nbf(this, &ShellData::assert_not_initialized),
    center(this, &ShellData::assert_not_initialized),
    atom_bfoff(this, &ShellData::assert_not_initialized),
    atom_shoff(this, &ShellData::assert_not_initialized),
    atom_nsh(this, &ShellData::assert_not_initialized),
    atom_nbf(this, &ShellData::assert_not_initialized),
    bfoff_in_atom(this, &ShellData::assert_not_initialized),
    shoff_in_atom(this, &ShellData::assert_not_initialized),
    atom_last_function(this, &ShellData::assert_not_initialized),
    atom_last_shell(this, &ShellData::assert_not_initialized),
    last_function(this, &ShellData::assert_not_initialized),
    atom_dfshoff(this, &ShellData::assert_not_initialized),
    atom_dfbfoff(this, &ShellData::assert_not_initialized),
    atom_dfnbf(this, &ShellData::assert_not_initialized),
    atom_dfnsh(this, &ShellData::assert_not_initialized),
    atom_df_last_function(this, &ShellData::assert_not_initialized),
    atom_df_last_shell(this, &ShellData::assert_not_initialized)
{

}

inline
ShellData::ShellData(
    int idx,
    Ref<GaussianBasisSet> basis,
    Ref<GaussianBasisSet> dfbasis
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
  assert(basis->nshell() > 0);
}

inline
ShellData::ShellData(
    const ShellData& other
) : index(other.index),
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
  assert(basis->nshell() > 0);
}

inline ShellData
shell_data(
    Ref<GaussianBasisSet> basis,
    int ish, Ref<GaussianBasisSet> dfbasis
)
{
  return ShellData(ish, basis, dfbasis);
}

inline ShellData&
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

//============================================================================//
// BasisFunctionData

inline
BasisFunctionData::BasisFunctionData()
  : index(-1),
    shell_index(this, &BasisFunctionData::assert_not_initialized),
    shell_nbf(this, &BasisFunctionData::assert_not_initialized),
    shell_bfoff(this, &BasisFunctionData::assert_not_initialized),
    center(this, &BasisFunctionData::assert_not_initialized),
    atom_bfoff(this, &BasisFunctionData::assert_not_initialized),
    atom_shoff(this, &BasisFunctionData::assert_not_initialized),
    atom_nsh(this, &BasisFunctionData::assert_not_initialized),
    atom_nbf(this, &BasisFunctionData::assert_not_initialized),
    bfoff_in_atom(this, &BasisFunctionData::assert_not_initialized),
    bfoff_in_shell(this, &BasisFunctionData::assert_not_initialized),
    shoff_in_atom(this, &BasisFunctionData::assert_not_initialized),
    atom_last_function(this, &BasisFunctionData::assert_not_initialized),
    atom_last_shell(this, &BasisFunctionData::assert_not_initialized),
    atom_dfshoff(this, &BasisFunctionData::assert_not_initialized),
    atom_dfbfoff(this, &BasisFunctionData::assert_not_initialized),
    atom_dfnbf(this, &BasisFunctionData::assert_not_initialized),
    atom_dfnsh(this, &BasisFunctionData::assert_not_initialized),
    atom_df_last_function(this, &BasisFunctionData::assert_not_initialized),
    atom_df_last_shell(this, &BasisFunctionData::assert_not_initialized)
{

}

inline
BasisFunctionData::BasisFunctionData(
    int idx,
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis
) : index(idx),
    basis(basis),
    dfbasis(dfbasis),
    shell_index(this, &BasisFunctionData::get_shell_index),
    shell_nbf(this, &BasisFunctionData::get_shell_nbf),
    shell_bfoff(this, &BasisFunctionData::get_shell_bfoff),
    center(this, &BasisFunctionData::get_center),
    atom_bfoff(this, &BasisFunctionData::get_atom_bfoff),
    atom_shoff(this, &BasisFunctionData::get_atom_shoff),
    atom_nsh(this, &BasisFunctionData::get_atom_nsh),
    atom_nbf(this, &BasisFunctionData::get_atom_nbf),
    bfoff_in_atom(this, &BasisFunctionData::get_bfoff_in_atom),
    bfoff_in_shell(this, &BasisFunctionData::get_bfoff_in_shell),
    shoff_in_atom(this, &BasisFunctionData::get_shoff_in_atom),
    atom_last_function(this, &BasisFunctionData::get_atom_last_function),
    atom_last_shell(this, &BasisFunctionData::get_atom_last_shell),
    atom_dfshoff(this, &BasisFunctionData::get_atom_dfshoff),
    atom_dfbfoff(this, &BasisFunctionData::get_atom_dfbfoff),
    atom_dfnbf(this, &BasisFunctionData::get_atom_dfnbf),
    atom_dfnsh(this, &BasisFunctionData::get_atom_dfnsh),
    atom_df_last_function(this, &BasisFunctionData::get_atom_df_last_function),
    atom_df_last_shell(this, &BasisFunctionData::get_atom_df_last_shell)
{

}

inline
BasisFunctionData::BasisFunctionData(
    const BasisFunctionData& other
) : index(other.index),
    basis(other.basis),
    dfbasis(other.dfbasis),
    shell_index(this, &BasisFunctionData::get_shell_index),
    shell_nbf(this, &BasisFunctionData::get_shell_nbf),
    shell_bfoff(this, &BasisFunctionData::get_shell_bfoff),
    center(this, &BasisFunctionData::get_center),
    atom_bfoff(this, &BasisFunctionData::get_atom_bfoff),
    atom_shoff(this, &BasisFunctionData::get_atom_shoff),
    atom_nsh(this, &BasisFunctionData::get_atom_nsh),
    atom_nbf(this, &BasisFunctionData::get_atom_nbf),
    bfoff_in_atom(this, &BasisFunctionData::get_bfoff_in_atom),
    bfoff_in_shell(this, &BasisFunctionData::get_bfoff_in_shell),
    shoff_in_atom(this, &BasisFunctionData::get_shoff_in_atom),
    atom_last_function(this, &BasisFunctionData::get_atom_last_function),
    atom_last_shell(this, &BasisFunctionData::get_atom_last_shell),
    atom_dfshoff(this, &BasisFunctionData::get_atom_dfshoff),
    atom_dfbfoff(this, &BasisFunctionData::get_atom_dfbfoff),
    atom_dfnbf(this, &BasisFunctionData::get_atom_dfnbf),
    atom_dfnsh(this, &BasisFunctionData::get_atom_dfnsh),
    atom_df_last_function(this, &BasisFunctionData::get_atom_df_last_function),
    atom_df_last_shell(this, &BasisFunctionData::get_atom_df_last_shell)
{

}

inline BasisFunctionData&
BasisFunctionData::operator=(const BasisFunctionData& other)
{
  index = other.index;
  basis = other.basis;
  dfbasis = other.dfbasis;
  shell_index.set_target_getter_setter(this, &BasisFunctionData::get_shell_index);
  shell_nbf.set_target_getter_setter(this, &BasisFunctionData::get_shell_nbf);
  shell_bfoff.set_target_getter_setter(this, &BasisFunctionData::get_shell_bfoff);
  center.set_target_getter_setter(this, &BasisFunctionData::get_center);
  atom_bfoff.set_target_getter_setter(this, &BasisFunctionData::get_atom_bfoff);
  atom_shoff.set_target_getter_setter(this, &BasisFunctionData::get_atom_shoff);
  atom_nsh.set_target_getter_setter(this, &BasisFunctionData::get_atom_nsh);
  atom_nbf.set_target_getter_setter(this, &BasisFunctionData::get_atom_nbf);
  bfoff_in_atom.set_target_getter_setter(this, &BasisFunctionData::get_bfoff_in_atom);
  bfoff_in_shell.set_target_getter_setter(this, &BasisFunctionData::get_bfoff_in_shell);
  shoff_in_atom.set_target_getter_setter(this, &BasisFunctionData::get_shoff_in_atom);
  atom_last_function.set_target_getter_setter(this, &BasisFunctionData::get_atom_last_function);
  atom_last_shell.set_target_getter_setter(this, &BasisFunctionData::get_atom_last_shell);
  atom_dfshoff.set_target_getter_setter(this, &BasisFunctionData::get_atom_dfshoff);
  atom_dfbfoff.set_target_getter_setter(this, &BasisFunctionData::get_atom_dfbfoff);
  atom_dfnbf.set_target_getter_setter(this, &BasisFunctionData::get_atom_dfnbf);
  atom_dfnsh.set_target_getter_setter(this, &BasisFunctionData::get_atom_dfnsh);
  atom_df_last_function.set_target_getter_setter(this, &BasisFunctionData::get_atom_df_last_function);
  atom_df_last_shell.set_target_getter_setter(this, &BasisFunctionData::get_atom_df_last_shell);
  return *this;
}


inline const BasisFunctionData
function_data(
    const Ref<GaussianBasisSet>& basis,
    int ish, const Ref<GaussianBasisSet>& dfbasis
)
{
  return BasisFunctionData(ish, basis, dfbasis);
}

//============================================================================//
// shell_iterator_wrapper

inline
shell_iter_wrapper::shell_iter_wrapper(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_shell,
    int last_shell
) : basis_(basis),
    dfbasis_(dfbasis),
    first_shell_(first_shell),
    last_shell_(last_shell)
{

}

inline ShellIter
shell_iter_wrapper::begin() const
{
  if(dfbasis_.nonnull()){
    ShellIter rv(basis_, dfbasis_, first_shell_);
    return rv;
  }
  else{
    ShellIter rv(basis_, first_shell_);
    return rv;
  }
}

inline ShellIter
shell_iter_wrapper::end() const
{
  if(dfbasis_.nonnull()){
    ShellIter rv(basis_, dfbasis_, last_shell_ == -1 ? basis_->nshell() : last_shell_ + 1);
    return rv;
  }
  else{
    ShellIter rv(basis_, last_shell_ == -1 ? basis_->nshell() : last_shell_ + 1);
    return rv;
  }
}



//============================================================================//

inline const shell_iter_wrapper
iter_shells_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(basis->shell_on_center(center, 0), basis, dfbasis);
  return shell_iter_wrapper(
      basis, dfbasis, sh.atom_shoff, sh.atom_last_shell
  );
}

inline const function_iter_wrapper
iter_functions_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(basis->shell_on_center(center, 0), basis, dfbasis);
  return function_iter_wrapper(
      basis, dfbasis, sh.atom_bfoff, sh.atom_last_function
  );
}

inline const function_iter_wrapper
iter_functions_in_shell(
    const Ref<GaussianBasisSet>& basis,
    int shell,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(shell, basis, dfbasis);
  return function_iter_wrapper(
      basis, dfbasis, sh.bfoff, sh.last_function
  );
}



} // end namespace sc



#endif /* _chemistry_qc_scf_cadf_iters_h */
