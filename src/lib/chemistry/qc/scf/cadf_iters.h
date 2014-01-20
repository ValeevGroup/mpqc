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
#include <boost/range.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <utility>
#include <type_traits>
#include <memory>
#include <iterator>

#define DEFAULT_TARGET_BLOCK_SIZE 100  // functions


namespace sc {

#define DUMP(expr) std::cout << #expr << " = " << (expr) << std::endl;
#define out_assert(a, op, b) assert(a op b || ((std::cout << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))

#define NOT_ASSIGNED -1

// Forward declarations
class ShellData;
class BasisFunctionData;
template<typename range_type> class ShellBlockIterator;
template<typename Iterator, typename ParentContainer> class IterableBasisElementData;

using int_range = decltype(boost::counting_range(1, 2));

struct ignored_argument { };

template <typename T>
class OptionalRefParameter {
    T* p_ = 0;

  public:

    OptionalRefParameter(T* p) : p_(p) { }

    OptionalRefParameter(const Ref<T>& rp) : p_(rp.nonnull() ? rp.pointer() : 0) { }

    operator T*() const { return p_; }
    operator const T*() const { return p_; }

};

typedef enum {
  NoRestrictions = 0,
  SameCenter = 1,
  SameAngularMomentum = 2
} BlockCompositionRequirement;

//############################################################################//

//############################################################################//

template <typename Iterator, typename DataContainer> struct BasisElementIteratorDereference;

class BasisElementData {

  public:

    BasisElementData(
        int index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis
    ) : index(index),
        basis(basis),
        dfbasis(dfbasis)
    { }

    BasisElementData() = delete;

    enum { NotAssigned = -1 };

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

    int index = NotAssigned;

    bool operator!=(const BasisElementData& other) const
    {
      return index != other.index;
    }

    virtual ~BasisElementData() { }

    int center = NotAssigned;
    int bfoff_in_atom = NotAssigned;

};

struct ShellData : public BasisElementData {

    template<typename Iterator> using with_iterator =
        BasisElementIteratorDereference<ShellData, Iterator>;

    ShellData(
        int index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0,
        int block_offset = NotAssigned
    ) : BasisElementData(index, basis, dfbasis),
        block_offset(NotAssigned)
    {
      init();
    }

    template <typename Iterator>
    ShellData(
        const with_iterator<Iterator>& deref
    ) : ShellData(deref.index, deref.basis, deref.dfbasis)
    { }

    ShellData() : BasisElementData(NotAssigned, 0, 0) { }


    int bfoff = NotAssigned;
    int nbf = NotAssigned;
    int atom_bfoff = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_nsh = NotAssigned;
    int atom_nbf = NotAssigned;
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

    // Reserved for future use
    int block_offset = NotAssigned;

    operator const int() const {
      out_assert(basis, !=, 0);
      out_assert(index, <, basis->nshell());
      out_assert(index, >=, 0);
      return index;
    }

    static int max_index(GaussianBasisSet* basis) {
      return basis->nshell();
    }

  protected:

    void init(){

      if(index == NotAssigned || index == basis->nshell()) return;

      out_assert(basis, !=, 0);
      out_assert(index, <, basis->nshell());
      out_assert(index, >=, 0);

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

    template<typename Iterator> using with_iterator =
        BasisElementIteratorDereference<BasisFunctionData, Iterator>;

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

    template <typename Iterator>
    BasisFunctionData(
        const with_iterator<Iterator>& deref
    ) : BasisFunctionData(deref.index, deref.basis, deref.dfbasis, deref.block_offset)
    { }

    BasisFunctionData()
      : BasisFunctionData(NotAssigned, 0, 0)
    { }

    int shell_index = NotAssigned;
    int shell_bfoff = NotAssigned;
    int bfoff_in_shell = NotAssigned;
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_bfoff = NotAssigned;
    int block_offset = NotAssigned;
    int bfoff_in_block = NotAssigned;

    operator const int() const { return index; }

    static int max_index(GaussianBasisSet* basis) {
      return basis->nbasis();
    }

  protected:

    void init()
    {
      if(index == NotAssigned || index == basis->nbasis()) return;

      out_assert(basis, !=, 0);
      out_assert(index, <, basis->nbasis());
      out_assert(index, >=, 0);

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


template <typename DataContainer, typename Iterator>
struct BasisElementIteratorDereference : DataContainer {

    const Iterator iterator;

  public:

    template<typename... Args>
    BasisElementIteratorDereference(
        Iterator iter,
        Args&&... args
    ) : DataContainer(std::forward<Args>(args)...), iterator(iter)
    { }

    //operator const Iterator(){ return iterator; }

};


//############################################################################//
//############################################################################//
//############################################################################//

enum { NotAssigned = -1 };

template <typename DataContainer, typename Iterator=int_range::iterator>
class basis_element_iterator
  : public boost::iterator_adaptor<
      basis_element_iterator<DataContainer, Iterator>,
      Iterator,
      boost::use_default,
      boost::bidirectional_traversal_tag,
      BasisElementIteratorDereference<DataContainer, Iterator>
    >
{

  public:

    typedef basis_element_iterator<DataContainer, Iterator> self_t;
    typedef boost::iterator_adaptor<
         self_t, Iterator, boost::use_default, boost::bidirectional_traversal_tag,
         BasisElementIteratorDereference<DataContainer, Iterator>
     > super_t;
    typedef BasisElementIteratorDereference<DataContainer, Iterator> value_reference_t;


  private:

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int block_offset = NotAssigned;

    friend class boost::iterator_core_access;

    value_reference_t dereference() const {
      auto base_spot = super_t::base();
      return value_reference_t(
        base_spot,
        int(*base_spot),
        basis, dfbasis,
        block_offset
      );
    }

  public:

    basis_element_iterator() : basis(0), dfbasis(0) { }

    basis_element_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        Iterator iter,
        int block_offset = NotAssigned
    ) : basis_element_iterator::iterator_adaptor_(iter),
        basis(basis),
        dfbasis(dfbasis),
        block_offset(block_offset)
    { }

};

template<typename Iterator=int_range::iterator>
    using shell_iterator = basis_element_iterator<ShellData, Iterator>;

template<typename Iterator=int_range::iterator>
    using function_iterator = basis_element_iterator<BasisFunctionData, Iterator>;


template<typename DataContainer, typename Iterator=int_range::iterator>
    using range_of = decltype(boost::make_iterator_range(
        basis_element_iterator<DataContainer, Iterator>(),
        basis_element_iterator<DataContainer, Iterator>()
    ));

enum {
  NoLastIndex = -1
};

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    OptionalRefParameter<GaussianBasisSet> dfbasis,
    int first_index,
    int last_index,
    int block_offset = NotAssigned
)
{
  const int actual_last_index = last_index == NoLastIndex ?
      DataContainer::max_index(basis) - 1 : last_index;
  auto base_range = boost::counting_range(first_index, actual_last_index + 1);
  return boost::make_iterator_range(
      basis_element_iterator<DataContainer>(basis, dfbasis, base_range.begin(), block_offset),
      basis_element_iterator<DataContainer>(basis, dfbasis, base_range.end(), block_offset)
  );

}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    int last_index
)
{
  return basis_element_range<DataContainer>(basis, 0, 0, last_index);
}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis = 0,
    int last_index = NoLastIndex
)
{
  return basis_element_range<DataContainer>(basis, dfbasis, 0, last_index);
}

template <typename DataContainer>
inline range_of<DataContainer>
basis_element_range(
    GaussianBasisSet* basis,
    int first_index,
    int last_index
)
{
  return basis_element_range<DataContainer>(basis, 0, first_index, last_index);
}


template <typename Iterator>
inline range_of<ShellData, Iterator>
shell_range(
    const ShellData::with_iterator<Iterator>& begin,
    const ShellData::with_iterator<Iterator>& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, end.iterator)
  );
}

template <typename Iterator, typename Iterator2>
inline range_of<ShellData, Iterator>
shell_range(
    const ShellData::with_iterator<Iterator>& begin,
    const Iterator2& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<ShellData, Iterator>(begin.basis, begin.dfbasis, Iterator(end))
  );
}

template <typename Iterator>
inline range_of<BasisFunctionData, Iterator>
function_range(
    const BasisFunctionData::with_iterator<Iterator>& begin,
    const BasisFunctionData::with_iterator<Iterator>& end
)
{
  return boost::make_iterator_range(
      basis_element_iterator<BasisFunctionData, Iterator>(begin.basis, begin.dfbasis, begin.iterator),
      basis_element_iterator<BasisFunctionData, Iterator>(begin.basis, begin.dfbasis, end.iterator)
  );
}

template <typename... Args>
inline range_of<ShellData>
shell_range(GaussianBasisSet* basis, Args&&... args) {
  return basis_element_range<ShellData>(basis, std::forward<Args>(args)...);
}

template <typename... Args>
inline range_of<BasisFunctionData>
function_range(GaussianBasisSet* basis, Args&&... args) {
  return basis_element_range<BasisFunctionData>(basis, std::forward<Args>(args)...);
}

inline range_of<BasisFunctionData>
function_range(const ShellData& ish) {
  return basis_element_range<BasisFunctionData>(
      ish.basis, ish.dfbasis,
      ish.bfoff, ish.last_function
  );
}


template<typename Iterator> class ShellBlockSkeleton;

template<typename Range = range_of<ShellData>>
class ShellBlockData {

    void init();

  public:

    enum { NotAssigned = -1 };

    Range shell_range;
    ShellData first_shell;
    ShellData last_shell;
    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int restrictions;

    ShellBlockData() : restrictions(-1) { }

    ShellBlockData(const ShellBlockSkeleton<Range>&);
    
    explicit ShellBlockData(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : ShellBlockData(
          sc::shell_range(basis, dfbasis, 0, basis->nshell() - 1),
          basis->nshell(), basis->nbasis(),
          NoRestrictions
        )
    { }

    ShellBlockData(
        Range sh_range,
        int nshell, int nbf,
        int requirements
    ) : shell_range(sh_range),
        first_shell(*shell_range.begin()),
        last_shell(shell_range.begin() == shell_range.end() ?
            *shell_range.end() : *(--shell_range.end())
        ),
        basis(first_shell.basis),
        dfbasis(first_shell.dfbasis),
        nshell(nshell), nbf(nbf),
        restrictions(requirements)
    { init(); }

    int nbf;
    int bfoff;
    int nshell;
    int last_function;
    int center = NotAssigned;
    int atom_bfoff = NotAssigned;
    int atom_shoff = NotAssigned;
    int atom_nsh = NotAssigned;
    int atom_nbf = NotAssigned;
    int bfoff_in_atom = NotAssigned;
    int shoff_in_atom = NotAssigned;
    int atom_last_function = NotAssigned;
    int atom_last_shell = NotAssigned;
    int atom_dfshoff = NotAssigned;
    int atom_dfbfoff = NotAssigned;
    int atom_dfnbf = NotAssigned;
    int atom_dfnsh = NotAssigned;
    int atom_df_last_function = NotAssigned;
    int atom_df_last_shell = NotAssigned;

    static int max_index(GaussianBasisSet* basis) { return basis->nshell(); }

};

template<typename Range = range_of<ShellData>>
class ShellBlockSkeleton {

  private:

    int nshell;
    int nbf;


    friend class ShellBlockData<Range>;

  public:

    Range shell_range;
    int first_index;
    int restrictions;

    ShellBlockSkeleton() : nshell(0), nbf(0), first_index(-1) { }

    ShellBlockSkeleton(const ShellBlockData<Range>& shbd)
      : shell_range(shbd.shell_range), nshell(shbd.nshell), nbf(shbd.nbf),
        restrictions(shbd.restrictions),
        first_index((*shell_range.begin()).index)
    { }

    ShellBlockSkeleton(Range range, int nsh, int nbas, int restrictions)
      : shell_range(range), nshell(nsh), nbf(nbas),
        restrictions(restrictions),
        first_index((*shell_range.begin()).index)
    { }


    template<typename OtherRange>
    bool operator==(const ShellBlockSkeleton<OtherRange>& other) const
    {
      return !(this->operator!=(other));
    }

    template<typename OtherRange>
    bool operator!=(const ShellBlockSkeleton<OtherRange>& other) const
    {
      return first_index != other.first_index
          or nshell != other.nshell
          or nbf != other.nbf;
    }

    template<typename OtherRange>
    bool operator<(const ShellBlockSkeleton<OtherRange>& other) const
    {
      const int my_index = (*shell_range.begin()).index;
      const int other_index = (*other.shell_range.begin()).index;
      if(my_index < other_index) return true;
      else if(my_index > other_index) return false;
      else if(nshell < other.nshell) return true;
      else if(nshell > other.nshell) return false;
      else if(nbf < other.nbf) return true;
      else return false;
    }

};

template<typename ShellIterator>
class shell_block_iterator
  : public boost::iterator_facade<
      shell_block_iterator<ShellIterator>,
      ShellBlockSkeleton<range_of<ShellData, ShellIterator>>,
      boost::forward_traversal_tag,
      ShellBlockData<range_of<ShellData, ShellIterator>>
    >
{
  public:

    typedef range_of<ShellData, ShellIterator> ShellRange;
    typedef shell_block_iterator<ShellIterator> self_type;
    typedef ShellBlockData<range_of<ShellData, ShellIterator>> value_reference;

  private:

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    int target_size;
    int restrictions;

    ShellRange all_shells;
    ShellBlockSkeleton<ShellRange> current_skeleton;

    void init();

    friend class boost::iterator_core_access;

    value_reference dereference() const {
      return value_reference(
          current_skeleton
      );
    }

    template <typename OtherIterator>
    bool equal(shell_block_iterator<OtherIterator> const& other) const
    {
      return current_skeleton == other.current_skeleton;
    }

    void increment()
    {
      all_shells = shell_range(*current_skeleton.shell_range.end(), *all_shells.end());
      init();
    }

  public:

    enum {
      NoLastShell = -1,
      NoMaximumBlockSize = -2
    };

    shell_block_iterator() : basis(0), dfbasis(0), restrictions(0), target_size(0) { }

    shell_block_iterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0,
        int first_index = 0,
        int last_index = NoLastShell,
        int requirements = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : basis(basis),
        dfbasis(dfbasis),
        restrictions(requirements),
        target_size(target_size),
        all_shells(shell_range(
            basis, dfbasis, first_index,
            last_index == NoLastShell ? ShellData::max_index(basis) - 1 : last_index
        ))
    {
      init();
    }

    static shell_block_iterator end_of_basis(GaussianBasisSet* basis) {
      return shell_block_iterator(basis, 0, basis->nshell(), basis->nshell()-1);
    }

    static shell_block_iterator end_with_last_index(int index, GaussianBasisSet* basis) {
      // The end iterator needs to start after the "last_index" that is passed in.  But
      //   because last_index gets incremented before being made into a counting range,
      //   the last_index parameter of the end iterator should be one less than the first.
      //   This is counterintuitive, but not meant to be exposed to the outside world anyway.
      return shell_block_iterator(basis, 0, index + 1, index);
    }

    const self_type& requiring(int restr) {
      restrictions = restr;
      init();
      return *this;
    }

};

// Type alias specialization for shell blocks
template<typename Iterator=int_range::iterator>
    using range_of_shell_blocks = decltype(boost::make_iterator_range(
        shell_block_iterator<Iterator>(),
        shell_block_iterator<Iterator>()
    ));

inline range_of_shell_blocks<>
shell_block_range(
    GaussianBasisSet* basis,
    OptionalRefParameter<GaussianBasisSet> dfbasis = 0,
    int first_index = 0,
    int last_index = NoLastIndex,
    int reqs = SameCenter,
    int target_size = DEFAULT_TARGET_BLOCK_SIZE
)
{
  typedef shell_block_iterator<int_range::iterator> result_t;
  const int actual_last_index = last_index == NoLastIndex ?
      basis->nshell() - 1 : last_index;
  return boost::make_iterator_range(
      result_t(basis, dfbasis, first_index,
          actual_last_index, reqs, target_size
      ),
      result_t::end_with_last_index(actual_last_index, basis)
  );

}

//############################################################################//

inline range_of<BasisFunctionData>
iter_functions_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  const int bfoff = basis->shell_to_function(shoff);
  return function_range(
      basis, dfbasis,
      bfoff, bfoff + basis->nbasis_on_center(center) - 1
  );
}

inline range_of<ShellData>
iter_shells_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_range(
      basis, dfbasis,
      shoff, shoff + basis->nshell_on_center(center) - 1
  );
}

inline range_of_shell_blocks<>
iter_shell_blocks_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0,
    int reqs = SameCenter,
    int target_size = DEFAULT_TARGET_BLOCK_SIZE
)
{
  const int shoff = basis->shell_on_center(center, 0);
  return shell_block_range(
      basis, dfbasis,
      shoff,
      shoff + basis->nshell_on_center(center) - 1,
      reqs|SameCenter,
      target_size
  );
}



} // end namespace sc

#include <chemistry/qc/scf/cadf_iters_impl.h>

#endif /* _chemistry_qc_scf_cadf_iters_h */
