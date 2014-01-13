//
// cadf_iters_impl.h
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

#ifndef _chemistry_qc_scf_cadf_iters_impl_h
#define _chemistry_qc_scf_cadf_iters_impl_h


#define out_assert(a, op, b) assert(a op b || ((std::cout << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))

namespace sc {

//============================================================================//
// ShellData

inline ShellData
shell_data(
    GaussianBasisSet* basis,
    int ish, GaussianBasisSet* dfbasis
)
{
  return ShellData(ish, basis, dfbasis);
}

inline const BasisFunctionData
function_data(
    const Ref<GaussianBasisSet>& basis,
    int ish, const Ref<GaussianBasisSet>& dfbasis
)
{
  return BasisFunctionData(ish, basis, dfbasis);
}

template <typename DataContainer>
inline typename detail::basis_iterator<DataContainer>::value_type
detail::basis_iterator<DataContainer>::begin() const
{
  return detail::basis_iterator<DataContainer>::value_type(first_index_, basis_, dfbasis_);
}

template <typename DataContainer>
inline typename detail::basis_iterator<DataContainer>::value_type
detail::basis_iterator<DataContainer>::end() const
{
  assert(last_index_ != NoLastIndex);
  return detail::basis_iterator<DataContainer>::value_type(last_index_, basis_, dfbasis_);
}

inline ShellData
shell_iterator::end() const
{
  return ShellData(last_index_ == NoLastIndex ? basis_->nshell() : last_index_ + 1, basis_, dfbasis_);
}

inline
function_iterator::function_iterator(const ShellData& ish)
  : detail::basis_iterator<BasisFunctionData>(
      ish.basis,
      ish.dfbasis,
      ish.bfoff,
      ish.last_function
    )
{ }

inline BasisFunctionData
function_iterator::begin() const
{
  return BasisFunctionData(first_index_, basis_, dfbasis_, block_offset);
}

inline BasisFunctionData
function_iterator::end() const
{
  return BasisFunctionData(last_index_ == -1 ? basis_->nbasis() : last_index_ + 1, basis_, dfbasis_);
}

//============================================================================//

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

inline void
ShellBlockData::init()
{
  // Treat the end iterator specially
  if(first_index == basis->nshell()){
    last_index = first_index;
    return;
  }
  //----------------------------------------//
  first_shell = ShellData(first_index, basis, dfbasis);
  nbf = 0;
  int nshell = 0;
  const int first_center = first_shell.center;
  auto& first_am = basis->shell(first_shell).am();
  for(auto ish : shell_iterator(basis, dfbasis, first_index, basis->nshell() - 1)){
    if(
        ((reqs & SameCenter) and ish.center != first_center) or
        ((reqs & SameAngularMomentum) and
            basis->shell(ish).am() != first_am)
    ){
      break;
    }
    else{
      ++nshell;
      nbf += ish.nbf;
      if(nbf >= target_size and target_size != -1) break;
    }
  }
  last_index = first_index + nshell - 1;
  last_shell = ShellData(last_index, basis, dfbasis);
  //----------------------------------------//
  bfoff = first_shell.bfoff;
  last_function = last_shell.last_function;
  if(reqs & SameCenter) {
    center = first_shell.center;
    atom_bfoff = first_shell.atom_bfoff;
    atom_shoff = first_shell.atom_shoff;
    atom_nsh = first_shell.atom_nsh;
    atom_nbf = first_shell.atom_nbf;
    bfoff_in_atom = first_shell.bfoff_in_atom;
    shoff_in_atom = first_shell.shoff_in_atom;
    atom_last_function = last_shell.atom_last_function;
    atom_last_shell = last_shell.atom_last_function;
    if(dfbasis != 0){
      atom_dfshoff = first_shell.atom_dfshoff;
      atom_dfbfoff = first_shell.atom_dfbfoff;
      atom_dfnbf = first_shell.atom_dfnbf;
      atom_dfnsh = first_shell.atom_dfnsh;
      atom_df_last_function = last_shell.atom_df_last_function;
      atom_df_last_shell = last_shell.atom_df_last_shell;
    }
  }
}

//============================================================================//

inline
shell_iterator::shell_iterator(
    const ShellBlockData& block
) : detail::basis_iterator<ShellData>(
      block.basis,
      block.dfbasis,
      block.first_shell.index,
      block.last_shell.index
    )
{ }

inline
function_iterator::function_iterator(
    const ShellBlockData& block
) : detail::basis_iterator<BasisFunctionData>(
      block.basis,
      block.dfbasis,
      block.bfoff,
      block.last_function
    ),
    block_offset(block.bfoff)
{ }

//============================================================================//

inline ShellBlockData
shell_block_iterator::begin() const
{
  return ShellBlockData(first_index_, basis_, dfbasis_, reqs, target_size);
}

inline ShellBlockData
shell_block_iterator::end() const
{
  const int end_shell = last_index_ == -1 ? basis_->nshell() : last_index_ + 1;
  return ShellBlockData(end_shell, basis_, dfbasis_, reqs, target_size);
}


//============================================================================//

} // end namespace sc


#endif /* _chemistry_qc_scf_cadf_iters_impl_h */
