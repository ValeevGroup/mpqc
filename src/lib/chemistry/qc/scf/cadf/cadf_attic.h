//
// cadf_attic.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jan 13, 2014
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

#ifndef CADF_ATTIC_H_
#define CADF_ATTIC_H_

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/logical/not.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/tuple/rem.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#define DEF_PROPERTY_use_properties(clname, type, condition, name, getter) \
      private:\
      inline const type BOOST_PP_CAT(get_, name)() const { assert(condition); return getter }\
      public:\
      property<clname, const type> name;
#define DEF_PROPERTY_no_use_properties(clname, type, condition, name, getter) type name = NotAssigned;
#define DEF_PROPERTY_impl(clname, type, condition, name, getter) \
  BOOST_PP_IF(BOOST_PP_CAT(clname, _USE_PROPERTIES), DEF_PROPERTY_use_properties, DEF_PROPERTY_no_use_properties)(clname, type, condition, name, getter)
#define DEF_PROPERTY(r,data,impl) \
    BOOST_PP_IF(BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 2), DEF_PROPERTY_impl, CONDITIONAL_PROPERTY_impl)(\
        BOOST_PP_SEQ_ELEM(0, data), BOOST_PP_SEQ_ELEM(1, data), BOOST_PP_SEQ_ELEM(2, data),\
        BOOST_PP_SEQ_ELEM(BOOST_PP_ADD(0, BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 3)), impl),\
        BOOST_PP_SEQ_ELEM(BOOST_PP_ADD(1, BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 3)), impl)\
    )
#define DEF_PROPERTY2(z,i,dataprops) DEF_PROPERTY2_impl(BOOST_PP_SEQ_ELEM(0,dataprops), BOOST_PP_SEQ_ELEM(i,BOOST_PP_SEQ_ELEM(1,dataprops)))
#define DEF_PROPERTY2_impl(data, impl) BOOST_PP_IF(BOOST_PP_CAT(BOOST_PP_SEQ_ELEM(0, data), _USE_PROPERTIES), DEF_PROPERTY_use_properties, DEF_PROPERTY_no_use_properties)(\
        BOOST_PP_SEQ_ELEM(0, data), BOOST_PP_SEQ_ELEM(1, data), BOOST_PP_SEQ_ELEM(2, data),\
        BOOST_PP_SEQ_ELEM(0, impl), BOOST_PP_SEQ_ELEM(1, impl)\
    )
#define CONDITIONAL_PROPERTY_impl(clname, type, old_conds, condition, props) BOOST_PP_REPEAT(BOOST_PP_SEQ_SIZE(props), DEF_PROPERTY2, ((clname)(type)((old_conds) && (condition)))(props))
#define CONDITIONAL_PROPERTIES(condition, props) (1)(condition)(props)

#define ASSERT_INITIALIZED(r,data,prop) assert(BOOST_PP_SEQ_ELEM(0,prop) != NotAssigned);
#define PROP_INIT_impl(check_these, pname, getter) BOOST_PP_SEQ_FOR_EACH(ASSERT_INITIALIZED, ~, check_these) pname = getter
#define PROP_INIT_cond_impl(check_these, cond, props) if(cond) { BOOST_PP_SEQ_FOR_EACH(PROP_INIT2, ~, props) }

#define PROP_INIT2(r,data,impl) BOOST_PP_SEQ_ELEM(0, impl) = BOOST_PP_SEQ_ELEM(1, impl)

#define PROP_INIT(r,data,i,impl) BOOST_PP_IF(BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 2), PROP_INIT_impl, PROP_INIT_cond_impl)(\
        BOOST_PP_SEQ_FIRST_N(i,data), \
        BOOST_PP_SEQ_ELEM(BOOST_PP_ADD(0, BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 3)), impl),\
        BOOST_PP_SEQ_ELEM(BOOST_PP_ADD(1, BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(impl), 3)), impl)\
)

#define PROPS_ASSIGN_nocond(clname, prop) BOOST_PP_SEQ_ELEM(0,prop).set_target_getter_setter(this, BOOST_PP_CAT(&clname::get_, BOOST_PP_SEQ_ELEM(0,prop)));
#define PROPS_ASSIGN_cond_impl(r, i, prop) BOOST_PP_SEQ_ELEM(0,prop).set_target(this);
#define PROPS_ASSIGN_cond(clname, prop) BOOST_PP_SEQ_FOR_EACH(PROPS_ASSIGN_cond_impl, ~, BOOST_PP_SEQ_ELEM(2,prop))
#define PROPS_ASSIGN_impl(r,data,i,prop) BOOST_PP_IF(BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(prop), 2), PROPS_ASSIGN_nocond, PROPS_ASSIGN_cond)(data, prop)
#define PROPS_ASSIGN(clname, props) BOOST_PP_SEQ_FOR_EACH_I(PROPS_ASSIGN_impl, clname, props)
#define NO_PROPS_ASSIGN(clname, props) init();
#define DEFINE_PROPERTIES(clname, type, props, init_function_name) \
    /* Create the member variables */ \
    BOOST_PP_SEQ_FOR_EACH_R(1, DEF_PROPERTY, (clname)(type)(true),  props) \
    /* create the init() function */ \
    BOOST_PP_IF(BOOST_PP_CAT(clname, _USE_PROPERTIES), \
        private: \
            void init_function_name() { \
              BOOST_PP_SEQ_FOR_EACH_I(PROPS_ASSIGN_impl, clname, props) \
            } \
        public: \
        , \
        private: \
            void init_function_name() { \
              BOOST_PP_SEQ_FOR_EACH_I( PROP_INIT, props, props) \
            } \
        public: \
    )

#define ShellData_USE_PROPERTIES 1
struct ShellData : public BasisElementData {


    DEFINE_PROPERTIES(ShellData, int,
        ((bfoff)(
            basis->shell_to_function(index);
        ))
        ((nbf)(
            basis->shell(index).nfunction();
        ))
        ((center)(
            basis->shell_to_center(index);
        ))
        ((atom_shoff)(
            basis->shell_on_center(center, 0);
        ))
        ((atom_bfoff)(
            basis->shell_to_function(atom_shoff);
        ))
        ((atom_nsh)(
            basis->nshell_on_center(center);
        ))
        ((atom_nbf)(
            basis->nbasis_on_center(center);
        ))
        ((shoff_in_atom)(
            index - atom_shoff;
        ))
        ((bfoff_in_atom)(
            bfoff - atom_bfoff;
        ))
        ((atom_last_function)(
            atom_bfoff + atom_nbf - 1;
        ))
        ((atom_last_shell)(
            atom_shoff + atom_nsh - 1;
        ))
        ((last_function)(
            bfoff + nbf - 1;
        ))
        (CONDITIONAL_PROPERTIES(dfbasis != 0,
            ((atom_dfshoff)(
                dfbasis->shell_on_center(center, 0);
            ))
            ((atom_dfbfoff)(
                dfbasis->shell_to_function(atom_dfshoff);
            ))
            ((atom_dfnsh)(
                dfbasis->nshell_on_center(center);
            ))
            ((atom_dfnbf)(
                dfbasis->nbasis_on_center(center);
            ))
            ((atom_df_last_function)(
                atom_dfbfoff + atom_dfnbf - 1; ))
            ((atom_df_last_shell)(
                atom_dfshoff + atom_dfnsh - 1;
            ))
        )),
        init
    );

    ShellData()
      : BasisElementData(NotAssigned, 0, 0)
    { }

    ShellData(
        int idx,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : BasisElementData(idx, basis, dfbasis)
    {
      init();
    }

    ShellData(
        const ShellData& other
    ) : BasisElementData(other.index, other.basis, other.dfbasis)
    {
      init();
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
      if(&other != this){
        init();
      }
      return *this;
    }

    operator int() { ASSERT_SHELL_BOUNDS; return index; }
    operator const int() const { ASSERT_SHELL_BOUNDS; return index; }

};

//============================================================================//

namespace detail {

template <typename DataContainer, typename Iterable>
class basis_iterable {

  protected:

    enum { NoLastIndex = -1, NoFirstIndex = -2 };

    GaussianBasisSet* basis_;
    GaussianBasisSet* dfbasis_;
    const Iterable& indices_;

  private:

    std::shared_ptr<Iterable> created_index_iterator_;

  public:

    typedef basis_iterable<DataContainer, Iterable> self_type;
    typedef decltype(Iterable().begin()) index_iterator_type;
    typedef IterableBasisElementData<index_iterator_type, DataContainer> value_type;

    //basis_iterable() : basis_(0), dfbasis_(0), indices_() { }

    basis_iterable(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis, //OptionalRefParameter<GaussianBasisSet> dfbasis,
        int first_index,
        int last_index
    ) : basis_(basis), dfbasis_(dfbasis),
        created_index_iterator_(std::make_shared<Iterable>(
            boost::counting_range(
              first_index, last_index == NoLastIndex ? DataContainer::max_index(basis_) : last_index
            )
        )),
        indices_(*created_index_iterator_)
    { }

    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
    ) : basis_iterable(basis, dfbasis, 0, NoLastIndex)
    { }

    template<typename iterable_like, typename=typename std::enable_if<
        std::is_base_of<Iterable, typename std::decay<iterable_like>::type>::value
      >::type
    >
    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        const OptionalRefParameter<GaussianBasisSet>& dfbasis,
        const iterable_like& indices
    ) : basis_(basis), dfbasis_(dfbasis), indices_(indices)
    { }

    template<typename iterable_like, typename=typename std::enable_if<
        std::is_base_of<Iterable, typename std::decay<iterable_like>::type>::value
      >::type
    >
    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        const iterable_like& indices,
        const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
    ) : basis_(basis), dfbasis_(dfbasis), indices_(indices)
    { }

    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        int last_index,
        const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
    ) : basis_iterable(basis, dfbasis, 0, last_index)
    { }

    template<typename int_like, typename=typename std::enable_if<
        std::is_convertible<int_like, int>::value
      >::type
    >
    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis,
        int_like last_index
    ) : basis_iterable(basis, dfbasis, 0, (int)last_index)
    { }

    basis_iterable(
        const Ref<GaussianBasisSet>& basis,
        int first_index,
        int last_index,
        const OptionalRefParameter<GaussianBasisSet>& dfbasis = 0
    ) : basis_iterable(basis, dfbasis, first_index, last_index)
    { }

    basis_iterable(
        const DataContainer& first,
        const DataContainer& last
    ) : basis_iterable(first.basis, first.dfbasis, Iterable(first, last))
    { }


    value_type begin() const;
    value_type last() const;
    value_type end() const;

};

} // end namespace detail

//============================================================================//

template <typename Iterable>
class shell_iterable : public detail::basis_iterable<ShellData, Iterable> {

  public:

    typedef shell_iterable<Iterable> self_type;
    typedef typename detail::basis_iterable<ShellData, Iterable> base_type;
    typedef decltype(Iterable().begin()) index_iterator_type;
    typedef IterableBasisElementData<index_iterator_type, ShellData> value_type;

    // Not sure if I need these anymore or not...
    shell_iterable(self_type&) = default;
    shell_iterable(const self_type&) = default;

    explicit shell_iterable(const ShellBlockIterator<Iterable>& block);

    using detail::basis_iterable<ShellData, Iterable>::basis_iterable;

};

/**
 * Requirements on ShellIncrementable:  must be an incrementable with basis and dfbasis members.
 * Should also not be convertible to a BasisFunctionData object.
 */
template<typename ShellIncrementable>
shell_iterable<boost::iterator_range<boost::counting_iterator<ShellIncrementable>>>
make_shell_iterable(ShellIncrementable first, ShellIncrementable last) {

  static_assert(not std::is_convertible<ShellIncrementable, BasisFunctionData>::value,
      "Can't construct a shell_iterable using BasisFunctionData objects."
  );

  return shell_iterable<boost::iterator_range<boost::counting_iterator<ShellIncrementable>>>(
      first.basis,
      first.dfbasis,
      boost::iterator_range<boost::counting_iterator<ShellIncrementable>>(
          first, last
      )
  );
}

typedef shell_iterable<int_range> shell_range;

//============================================================================//

template <typename Iterable = int_range>
class function_iterable : public detail::basis_iterable<BasisFunctionData, Iterable> {

    static constexpr int NotAssigned = -1;

    int block_offset = NotAssigned;

  public:

    typedef function_iterable<Iterable> self_type;
    typedef typename detail::basis_iterable<BasisFunctionData, Iterable> base_type;
    typedef decltype(Iterable().begin()) index_iterator_type;
    typedef IterableBasisElementData<index_iterator_type, BasisFunctionData> value_type;

    // Not sure if I need these anymore or not...
    function_iterable() = default;
    function_iterable(function_iterable&) = default;
    function_iterable(const function_iterable&) = default;

    explicit function_iterable(const ShellData&);

    explicit function_iterable(const ShellBlockIterator<Iterable>&);

    using detail::basis_iterable<BasisFunctionData, Iterable>::basis_iterable;

};

typedef function_iterable<int_range> function_range;

//============================================================================//

template<typename Iterable>
class shell_block_iterable : public detail::basis_iterable<ShellBlockIterator<Iterable>, Iterable> {

    int target_size = DEFAULT_TARGET_BLOCK_SIZE;

    // Composed using bitwise or of BlockCompositionRequirement enums
    int reqs = SameCenter;

  public:

    typedef shell_block_iterable<Iterable> self_type;
    typedef typename detail::basis_iterable<ShellBlockIterator<Iterable>, Iterable> base_type;
    typedef decltype(Iterable().begin()) index_iterator_type;
    typedef IterableBasisElementData<index_iterator_type, ShellBlockIterator<Iterable>> value_type;

    shell_block_iterable() = default;
    shell_block_iterable(self_type&) = default;
    shell_block_iterable(const self_type&) = default;

    using detail::basis_iterable<ShellBlockIterator<Iterable>, Iterable>::basis_iterable;

    const shell_block_iterable& requiring(int in_reqs) {
      reqs = in_reqs;
      return *this;
    }

    const shell_block_iterable& with_target_size(int new_target_size) {
      target_size = new_target_size;
      return *this;
    }

    value_type begin() const;
    value_type last() const;
    value_type end() const;

};

typedef shell_block_iterable<int_range> shell_block_range;

template<typename Iterator, typename ParentContainer>
class IterableBasisElementData : public ParentContainer {

    Iterator spot;

  public:

    typedef IterableBasisElementData<Iterator, ParentContainer> self_type;
    typedef boost::forward_traversal_tag iterator_category;
    typedef typename Iterator::difference_type difference_type;

    template<typename... Args>
    IterableBasisElementData(
        Iterator index_iter,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        Args&&... args
    ) : spot(index_iter), ParentContainer(*index_iter, basis, dfbasis, std::forward<Args>(args)...)
    { }

    const self_type& operator++()
    {
      ++spot;
      ParentContainer::index = *spot;
      ParentContainer::init();
      return *this;
    }

    const self_type& operator*() const {
      return *this;
    }

};

template<typename range_type>
class ShellBlockIterator {

  public:

    static constexpr int NoMaximumBlockSize = -20;

    // Forward declaration of inner class
    class Skeleton;

    // Construct a ShellBlockIterator that spans an entire basis
    ShellBlockIterator(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    ) : BasisElementData(NotAssigned, basis, dfbasis),
        possible_shell_iter(basis, dfbasis, first_shell, basis->nshell()),
        target_size(NoMaximumBlockSize),
        reqs(NoRestrictions),
        first_shell(possible_shell_iter.begin()),
        last_shell(possible_shell_iter.end()), // for now; init will fix this
        shell_iter(first_shell, last_shell) // for now; init() will fix this
    {
      init();
    }

    ShellBlockIterator(
        int first_index,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int reqs = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : BasisElementData(NotAssigned, basis, dfbasis),
        possible_shell_iter(basis, dfbasis, first_index, basis->nshell()),
        target_size(target_size),
        reqs(reqs),
        first_shell(possible_shell_iter.begin()),
        last_shell(possible_shell_iter.end()), // for now; init will fix this
        shell_iter(first_shell, last_shell) // for now; init() will fix this
    {
      init();
    }

    template<typename iterable_like, typename=typename std::enable_if<
        std::is_base_of<Iterable, typename std::decay<iterable_like>::type>::value
        && std::is_convertible<iterable_like, Iterable>::value
      >::type>
    ShellBlockIterator(
        const iterable_like& iterable,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int reqs = SameCenter,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ) : BasisElementData(NotAssigned, basis, dfbasis),
        possible_shell_iter(basis, dfbasis, iterable),
        target_size(target_size),
        reqs(reqs),
        first_shell(possible_shell_iter.begin()),
        last_shell(possible_shell_iter.end()), // for now; init() will fix this
        shell_iter(first_shell, last_shell) // for now; init() will fix this
    {
      init();
    }


    ShellBlockIterator(const typename ShellBlockIterator<Iterable>::Skeleton&);

    /**
     * Note: not a full != operator; just for the purposes of range iteration!
     *   Two blocks with the same first index and the same size will compare
     *   equal, regardless of whether they contain the same indices.
     */
    bool operator!=(const ShellBlockIterator& other) const
    {
      return first_shell.index != other.first_shell.index or nshell != other.nshell;
    }

    const ShellBlockIterator& operator*() const { return *this; }

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;
    shell_iterable<Iterable> shell_iter;
    typename shell_iterable<Iterable>::value_type first_shell;
    typename shell_iterable<Iterable>::value_type last_shell;

    int nbf;
    int bfoff;
    int nshell;
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

    static int max_index(GaussianBasisSet* basis) { return basis->nshell(); }

  protected:

    shell_iterable<Iterable> possible_shell_iter;
    int first_index;
    int target_size;
    int reqs;

    void init();

  public:

    class Skeleton {

      public:

        int first_index;
        int size;
        int reqs;
        GaussianBasisSet* basis;
        GaussianBasisSet* dfbasis;

        Skeleton(
            const ShellBlockIterator& data
        ) : first_index(data.first_index),
            size(data.nbf),
            basis(data.basis),
            dfbasis(data.dfbasis),
            reqs(data.reqs)
        { }

        bool operator<(const ShellBlockIterator::Skeleton& other) const {
          return first_index < other.first_index
              or (first_index == other.first_index and size < other.size);
        }

        bool operator!=(const ShellBlockIterator::Skeleton& other) const {
          return first_index != other.first_index or size != other.size;
        }

    };

};

typedef ShellBlockIterator<int_range> ShellBlockData;

template<typename Iterator, typename Iterable>
class IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>> : public ShellBlockIterator<Iterable> {

    typedef IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>> self_type;
    typedef ShellBlockIterator<Iterable> parent_type;

    using ShellBlockIterator<Iterable>::basis;
    using ShellBlockIterator<Iterable>::dfbasis;
    using ShellBlockIterator<Iterable>::nshell;
    using ShellBlockIterator<Iterable>::possible_shell_iter;
    using ShellBlockIterator<Iterable>::nbf;
    using ShellBlockIterator<Iterable>::target_size;
    using ShellBlockIterator<Iterable>::last_shell;
    using ShellBlockIterator<Iterable>::first_shell;
    using ShellBlockIterator<Iterable>::reqs;
    using ShellBlockIterator<Iterable>::NoMaximumBlockSize;
    using ShellBlockIterator<Iterable>::init;

  protected:

    Iterator spot;

  public:

    // TODO avoid nesting of boost::counting_range types?
    template<typename... Args>
    IterableBasisElementData(
        const Iterable& index_iter,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        Args&&... args
    ) : spot(index_iter.begin()),
        ShellBlockIterator<Iterable>(
            index_iter,
            basis, dfbasis,
            std::forward<Args>(args)...
        )
    { }

    const self_type& operator++();

    const ShellBlockIterator<Iterable>& operator*() const { return *this; }

    const self_type& advance_to_last_block();

};

template <typename DataContainer, typename Iterable>
inline typename detail::basis_iterable<DataContainer, Iterable>::value_type
detail::basis_iterable<DataContainer, Iterable>::begin() const
{
  return value_type(indices_.begin(), basis_, dfbasis_);
}

template <typename DataContainer, typename Iterable>
inline typename detail::basis_iterable<DataContainer, Iterable>::value_type
detail::basis_iterable<DataContainer, Iterable>::last() const
{
  auto last = indices_.end();
  --last;
  return value_type(last, basis_, dfbasis_);
}

template <typename DataContainer, typename Iterable>
inline typename detail::basis_iterable<DataContainer, Iterable>::value_type
detail::basis_iterable<DataContainer, Iterable>::end() const
{
  return value_type(indices_.end(), basis_, dfbasis_);
}

template <typename Iterable>
inline typename shell_block_iterable<Iterable>::value_type
shell_block_iterable<Iterable>::begin() const
{
  return value_type(
      base_type::indices_,
      base_type::basis_, base_type::dfbasis_,
      reqs, target_size
   );
}

template <typename Iterable>
inline typename shell_block_iterable<Iterable>::value_type
shell_block_iterable<Iterable>::end() const
{
  return value_type(
      Iterable(base_type::indices_.end(), base_type::indices_.end()),
      base_type::basis_, base_type::dfbasis_,
      reqs, target_size
  );
}

template <typename Iterable>
inline typename shell_block_iterable<Iterable>::value_type
shell_block_iterable<Iterable>::last() const
{
  value_type rv(
      base_type::indices_.begin(),
      base_type::indices_.end(),
      base_type::basis_, base_type::dfbasis_,
      reqs, target_size
  );
  rv.advance_to_last_block();
  return rv;
}

template<typename Iterator, typename Iterable>
inline const IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>>&
IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>>::advance_to_last_block()
{

  //----------------------------------------//
  // Treat the end iterator specially
  if(possible_shell_iter.begin() == possible_shell_iter.end()){
    throw ProgrammingError(
        "Already at end; can't advance to last block, which is before end.",
        __FILE__, __LINE__
    );
  }
  auto iter_tmp = possible_shell_iter;
  const auto& last_possible_shell = possible_shell_iter.last();
  //----------------------------------------//
  while(iter_tmp.begin() != iter_tmp.end()) {
    first_shell = iter_tmp.begin();
    nbf = 0;
    nshell = 0;
    //----------------------------------------//
    const int first_center = first_shell.center;
    auto& first_am = basis->shell(first_shell).am();
    for(auto ish : possible_shell_iter){
      if(
          // Same center condition
          ((reqs & SameCenter) and ish.center != first_center)
          or
          // Same angular momentum condition
          ((reqs & SameAngularMomentum) and
              basis->shell(ish).am() != first_am) or
          // Maximum block size overflow condition
          (target_size != NoMaximumBlockSize and nbf >= target_size)
      ){
        // Should always have at least one shell, unless I mess up
        //   the conditional expression above.  Assert this fact:
        assert(nshell);
        // We're done making the block.  `ish` belongs to the next block.
        break;
      }
      else{
        // This could be the last shell.  Overwrite the last shell member
        //   everytime we get here
        last_shell = ish;
        ++nshell;
        nbf += ish.nbf;
      }
    }
    auto new_first_shell = last_shell;
    new_first_shell++;
    iter_tmp = make_shell_iterable(new_first_shell, last_possible_shell);
  }
  //----------------------------------------//
  possible_shell_iter = make_shell_iterable(first_shell, last_possible_shell);
  init();
}

template<typename Iterable>
inline
ShellBlockIterator<Iterable>::ShellBlockIterator(
    const typename ShellBlockIterator<Iterable>::Skeleton& sk
) : ShellBlockIterator<Iterable>(sk.first_index, sk.basis, sk.dfbasis, sk.reqs, sk.size)
{ }


template<typename Iterator, typename Iterable>
inline const IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>>&
IterableBasisElementData<Iterator, ShellBlockIterator<Iterable>>::operator++()
{
  spot += nshell;
  possible_shell_iter = make_shell_iter(spot, possible_shell_iter.last());
  init();
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// compute K attic

// three body parts
// Compute B
              /* Fail-safe explicit loops version
              //for(auto sigma : function_range(obs)) {
              //  B_mus[mu.bfoff_in_shell].row(sigma) += 2.0 *
              //      D.row(sigma).segment(jsh.bfoff, jsh.nbf) *
              //      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              //  //for(auto X : function_range(Xblk)) {
              //  //  B_mus[mu.bfoff_in_shell](sigma, X.bfoff_in_block) += 2.0 *
              //  //      D.row(sigma).segment(jsh.bfoff, jsh.nbf) *
              //  //      g3.col(X.bfoff_in_block).segment(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              //  //  //for(auto rho : function_range(jsh)) {
              //  //  //  B_mus[mu.bfoff_in_shell](sigma, X.bfoff_in_block) += 2.0 *
              //  //  //      g3(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell, X.bfoff_in_block)
              //  //  //      * D(rho, sigma);
              //  //  //}
              //  //
              //  //}
              //}
               */
// K Contributeions

            /* Failsafe version:
            //for(auto nu : iter_functions_on_center(obs, Y.center)) {
            //  //for(auto sigma : function_range(obs)) {
            //  //  Kt_part(mu, nu) -= 0.5 * C_Y(nu.bfoff_in_atom, sigma) * B_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
            //  //}
            //  Kt_part(mu, nu) -= 0.5 * C_Y.row(nu.bfoff_in_atom) * B_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell);
            //}
            //----------------------------------------//
            //for(auto nu : function_range(obs)) {
            //  if(nu.center != Y.center) {
            //    //for(auto sigma : iter_functions_on_center(obs, Y.center)) {
            //    //  Kt_part(mu, nu) -= 0.5 * C_Y(sigma.bfoff_in_atom, nu) * B_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
            //    //}
            //    Kt_part(mu, nu) -= 0.5 * C_Y.col(nu).transpose() * B_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell).segment(
            //      obs->shell_to_function(obs->shell_on_center(Y.center, 0)),
            //      obs->nbasis_on_center(Y.center)
            //    );
            //  }
            //}
              */
              //----------------------------------------//
          #if OLD_K
          for(auto ksh : shell_range(obs, dfbs_)) {
            for(auto lsh : iter_significant_partners(ksh)) {
              if(ksh.center != Xblk.center && lsh.center != Xblk.center) continue;
              for(auto nu : function_range(ksh)) {
                for(auto sigma : function_range(lsh)) {
                  CoefContainer c_ptr;
                  if(lsh <= ksh) {
                    IntPair nu_sigma(nu, sigma);
                    assert(coefs_.find(nu_sigma) != coefs_.end());
                    auto coef_pair = coefs_[nu_sigma];
                    c_ptr = ksh.center == Xblk.center ? coef_pair.first : coef_pair.second;
                  }
                  else {
                    IntPair nu_sigma(sigma, nu);
                    assert(coefs_.find(nu_sigma) != coefs_.end());
                    auto coef_pair = coefs_[nu_sigma];
                    c_ptr = lsh.center == Xblk.center ? coef_pair.first : coef_pair.second;
                  }
                  auto& C = *c_ptr;
                  for(auto X : function_range(Xblk)) {
                    for(auto mu : function_range(ish)) {
                      Kt_part(mu, nu) += C[X.bfoff_in_atom] * B_mus[mu.bfoff_in_shell](sigma, X.bfoff_in_block);
                    } // end loop over mu in ish
                  } // end loop over X in Xblk
                } // end loop over sigma in lsh
              } // end loop over nu in ksh
            } // end loop over lsh
          } // end loop over ksh
          #endif

// Full OLD_K section
        #if OLD_K
        ShellData jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          const double pf = 2.0;
          /*-----------------------------------------------------*/
          /* Setup                                          {{{2 */ #if 2 // begin fold
          std::vector<Eigen::MatrixXd> A_mus(ish.nbf);
          std::vector<Eigen::MatrixXd> A_rhos((ish == jsh) ? 0 : (int)jsh.nbf);
          std::vector<Eigen::MatrixXd> B_mus(ish.nbf);
          std::vector<Eigen::MatrixXd> B_rhos((ish == jsh) ? 0 : (int)jsh.nbf);
          std::vector<Eigen::MatrixXd> dt_mus(ish.nbf);
          std::vector<Eigen::MatrixXd> dt_rhos((ish == jsh) ? 0 : (int)jsh.nbf);
          std::vector<Eigen::MatrixXd> gt_mus(ish.nbf);
          std::vector<Eigen::MatrixXd> gt_rhos((ish == jsh) ? 0 : (int)jsh.nbf);
          for(auto mu : function_range(ish)){
            A_mus[mu.bfoff_in_shell].resize(nbf, dfnbf);
            A_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
            B_mus[mu.bfoff_in_shell].resize(nbf, dfnbf);
            B_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
            dt_mus[mu.bfoff_in_shell].resize(nbf, dfnbf);
            dt_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
            gt_mus[mu.bfoff_in_shell].resize(nbf, dfnbf);
            gt_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
          }
          if(ish != jsh){
            for(auto rho : function_range(jsh)){
              A_rhos[rho.bfoff_in_shell].resize(nbf, dfnbf);
              A_rhos[rho.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
              B_rhos[rho.bfoff_in_shell].resize(nbf, dfnbf);
              B_rhos[rho.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
              dt_rhos[rho.bfoff_in_shell].resize(nbf, dfnbf);
              dt_rhos[rho.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
              gt_rhos[rho.bfoff_in_shell].resize(nbf, dfnbf);
              gt_rhos[rho.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, dfnbf);
            }
          }
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Compute B for functions in ish, jsh            {{{2 */ #if 2 // begin fold
          if(ithr == 0) timer.enter("build B");
          for(int kshdf = 0; kshdf < dfbs_->nshell(); ++kshdf){
            const int kshdf_bfoff = dfbs_->shell_to_function(kshdf);
            const int kshdf_nbf = dfbs_->shell(kshdf).nfunction();
            //----------------------------------------//
            // Get the integrals for ish, jsh, ksh
            if(ithr == 0) timer.enter("compute ints");
            std::shared_ptr<Eigen::MatrixXd> g_part = ints_to_eigen(
                ish, jsh, kshdf,
                eris_3c_[ithr],
                coulomb_oper_type_
            );
            //----------------------------------------//
            // B_mus
            if(ithr == 0) timer.change("B_mus");
            for(int ibf = 0; ibf < ish.nbf; ++ibf){
              B_mus[ibf].middleCols(kshdf_bfoff, kshdf_nbf) +=
                  pf * D.middleCols(jsh.bfoff, jsh.nbf)
                    * g_part->middleRows(ibf * jsh.nbf, jsh.nbf);
            } // end loop over mu
            //----------------------------------------//
            if(ithr == 0) timer.change("B_rhos");
            if(ish != jsh) {
              // Do the contribution to the other way around
              for(int ibf = 0, mu = ish.bfoff; ibf < ish.nbf; ++ibf, ++mu){
                for(int jbf = 0; jbf < jsh.nbf; ++jbf){
                  //----------------------------------------//
                  B_rhos[jbf].middleCols(kshdf_bfoff, kshdf_nbf) +=
                      pf * D.col(mu) * g_part->row(ibf * jsh.nbf + jbf);
                  //----------------------------------------//
                }
              }
            } // end if ish != jsh
            //----------------------------------------//
            if(ithr == 0) timer.exit();
          } // end loop over kshdf
          //----------------------------------------//
          if(xml_debug) {
            for(auto mu : function_range(ish)){
              write_as_xml("B_part", B_mus[mu.bfoff_in_shell], std::map<std::string, int>{ {"mu", mu}, {"jsh", jsh} } );
            }
            if(ish != jsh){
              for(auto rho : function_range(jsh)){
                write_as_xml("B_part", B_rhos[rho.bfoff_in_shell], std::map<std::string, int>{ {"mu", rho}, {"jsh", jsh} } );
              }
            }
          }
          //----------------------------------------//
          if(ithr == 0) timer.exit("build B");
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Compute d_tilde and g_tilde for bfs in ish,jsh {{{2 */ #if 2 // begin fold
          //if(ithr == 0) timer.enter("build d_tilde");
          //for(auto mu : function_range(ish)){
          //  //----------------------------------------//
          //  // Compute d_tilde for a given mu
          //  for(auto rho : function_range(jsh)){
          //    //----------------------------------------//
          //    // Get the coefficients
          //    IntPair mu_rho(mu, rho);
          //    auto cpair = coefs_[mu_rho];
          //    VectorMap& Ca = *(cpair.first);
          //    VectorMap& Cb = *(cpair.second);
          //    //----------------------------------------//
          //    // dt_mus
          //    // column of D times Ca as row vector
          //    dt_mus[mu.bfoff_in_shell].middleCols(ish.atom_dfbfoff, ish.atom_dfnbf) +=
          //        pf * D.row(rho).transpose() * Ca.transpose();
          //    if(ish.center != jsh.center) {
          //      dt_mus[mu.bfoff_in_shell].middleCols(jsh.atom_dfbfoff, jsh.atom_dfnbf) +=
          //          pf * D.row(rho).transpose() * Cb.transpose();
          //    }
          //    //----------------------------------------//
          //    if(ish != jsh) {
          //      // dt_rhos
          //      // same thing, utilizing C_{mu,rho} = C_{rho,mu} (accounting for ordering, etc.)
          //      dt_rhos[rho.bfoff_in_shell].middleCols(ish.atom_dfbfoff, ish.atom_dfnbf) +=
          //          pf * D.row(mu).transpose() * Ca.transpose();
          //      if(ish.center != jsh.center) {
          //        dt_rhos[rho.bfoff_in_shell].middleCols(jsh.atom_dfbfoff, jsh.atom_dfnbf) +=
          //            pf * D.row(mu).transpose() * Cb.transpose();
          //      }
          //    }
          //    //----------------------------------------//
          //  } // end loop over rho
          //  //----------------------------------------//
          //} // end loop over mu
          //if(ithr == 0) timer.exit("build d_tilde");
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Get g_tilde contribution for bfs in ish,jsh    {{{2 */ #if 2 // begin fold
          //if(ithr == 0) timer.enter("build g");
          /* TODO optimally, this needs to compute blocks of several shells at once and contract these blocks
          for(auto Ysh : shell_range(dfbs_)) {
            for(auto Xsh : shell_range(dfbs_, 0, Ysh)) {
              if(ithr == 0) timer.enter("compute ints");
              const double dfpf = Xsh == Ysh ? 1.0 : 2.0;
              Eigen::MatrixXd& g2_part = *ints_to_eigen(
                Xsh, Ysh,
                eris_2c_[ithr],
                coulomb_oper_type_
              );
              //----------------------------------------//
              if(ithr == 0) timer.change("gt_mu");
              for(auto mu : function_range(ish)) {
                gt_mus[mu.bfoff_in_shell].middleCols(Ysh.bfoff, Ysh.nbf).noalias() +=
                    dfpf * dt_mus[mu.bfoff_in_shell].middleCols(Xsh.bfoff, Xsh.nbf) * g2_part;
              }
              //----------------------------------------//
              if(ish != jsh){
                if(ithr == 0) timer.change("gt_rho");
                for(auto rho : function_range(jsh)) {
                  gt_rhos[rho.bfoff_in_shell].middleCols(Ysh.bfoff, Ysh.nbf).noalias() +=
                      dfpf * dt_rhos[rho.bfoff_in_shell].middleCols(Xsh.bfoff, Xsh.nbf) * g2_part;
                } // end loop over rho
              }
              //----------------------------------------//
              if(ithr == 0) timer.exit();
            } // end loop over Xsh
          } // end loop over Ysh
          */
          //if(ithr == 0) timer.enter("gt_mu");
          //for(int mu = 0; mu < ish.nbf; ++mu){
          //  //gt_mus[mu].noalias() += dt_mus[mu] * g2;
          //  gt_mus[mu] = dt_mus[mu] * g2;
          //} // end loop over mu
          ////----------------------------------------//
          //if(ish != jsh){
          //  if(ithr == 0) timer.change("gt_rho");
          //  for(int rho = 0; rho < jsh.nbf; ++rho) {
          //    //gt_rhos[rho].noalias() += dt_rhos[rho] * g2;
          //    gt_rhos[rho] = dt_rhos[rho] * g2;
          //  } // end loop over rho
          //}
          //----------------------------------------//
          // Subtract off the term three contribution to B_mus and B_rhos (result is called A in notes)
          //if(ithr == 0) timer.enter("compute A");
          //for(auto mu : function_range(ish)) {
          //  A_mus[mu.bfoff_in_shell] = B_mus[mu.bfoff_in_shell] - 0.5 * dt_mus[mu.bfoff_in_shell] * g2; //gt_mus[mu.bfoff_in_shell];
          //}
          //if(ish != jsh) {
          //  for(auto rho : function_range(jsh)) {
          //    A_rhos[rho.bfoff_in_shell] = B_rhos[rho.bfoff_in_shell] - 0.5 * dt_rhos[rho.bfoff_in_shell] * g2; // gt_rhos[rho.bfoff_in_shell];
          //  }
          //}
          //if(ithr == 0) timer.exit("compute A");
          //if(ithr == 0) timer.exit("build g");
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Compute K contributions                        {{{2 */ #if 2 // begin fold
          if(ithr == 0) timer.enter("K contributions");
          //if(ithr == 0) timer.enter("misc");
          for(auto Y : function_range(dfbs_)) {
            Eigen::MatrixXd& C_Y = coefs_transpose_[Y];
            const int obs_atom_bfoff = obs->shell_to_function(obs->shell_on_center(Y.center, 0));
            const int obs_atom_nbf = obs->nbasis_on_center(Y.center);
            for(auto mu : function_range(ish)) {
              // B_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
              // C_Y is (Y.{obs_}atom_nbf x nbf)
              // result should be (Y.{obs_}atom_nbf x 1)
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() +=
                  C_Y * B_mus[mu.bfoff_in_shell].col(Y);
              // The sigma <-> nu term
              Kt_part.row(mu).transpose() += C_Y.transpose()
                  * B_mus[mu.bfoff_in_shell].col(Y).segment(obs_atom_bfoff, obs_atom_nbf);
              // Add back in the nu.center == Y.center part
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                  C_Y.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                  * B_mus[mu.bfoff_in_shell].col(Y).segment(obs_atom_bfoff, obs_atom_nbf);
              //----------------------------------------//
              /* Failsafe version:
              //for(auto nu : iter_functions_on_center(obs, Y.center)) {
              //  //for(auto sigma : function_range(obs)) {
              //  //  Kt_part(mu, nu) -= 0.5 * C_Y(nu.bfoff_in_atom, sigma) * B_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
              //  //}
              //  Kt_part(mu, nu) -= 0.5 * C_Y.row(nu.bfoff_in_atom) * B_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell);
              //}
              //----------------------------------------//
              //for(auto nu : function_range(obs)) {
              //  if(nu.center != Y.center) {
              //    //for(auto sigma : iter_functions_on_center(obs, Y.center)) {
              //    //  Kt_part(mu, nu) -= 0.5 * C_Y(sigma.bfoff_in_atom, nu) * B_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
              //    //}
              //    Kt_part(mu, nu) -= 0.5 * C_Y.col(nu).transpose() * B_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell).segment(
              //      obs->shell_to_function(obs->shell_on_center(Y.center, 0)),
              //      obs->nbasis_on_center(Y.center)
              //    );
              //  }
              //}
              */
              //----------------------------------------//
            }
            if(ish != jsh){
              for(auto rho : function_range(jsh)) {
                // B_rhos[rho.bfoff_in_shell] is (nbf x Ysh.nbf)
                // C_Y is (Y.{obs_}atom_nbf x nbf)
                // result should be (Y.{obs_}atom_nbf x 1)
                Kt_part.row(rho).segment(obs_atom_bfoff, obs_atom_nbf).transpose() +=
                    C_Y * B_rhos[rho.bfoff_in_shell].col(Y);
                // The sigma <-> nu term
                Kt_part.row(rho).transpose() += C_Y.transpose()
                    * B_rhos[rho.bfoff_in_shell].col(Y).segment(obs_atom_bfoff, obs_atom_nbf);
                // Add back in the nu.center == Y.center part
                Kt_part.row(rho).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                    C_Y.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                    * B_rhos[rho.bfoff_in_shell].col(Y).segment(obs_atom_bfoff, obs_atom_nbf);
                //----------------------------------------//
                /* Failsafe version:
                //for(auto nu : function_range(obs)) {
                //  for(auto sigma : function_range(obs)) {
                //    if(nu.center == Y.center){
                //      Kt_part(rho, nu) -= 0.5 * C_Y(nu.bfoff_in_atom, sigma) * dt_rhos[rho.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                //    }
                //    else if(sigma.center == Y.center){
                //      Kt_part(rho, nu) -= 0.5 * C_Y(sigma.bfoff_in_atom, nu) * dt_rhos[rho.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                //    }
                //  }
                //}
                */
                //----------------------------------------//
              }
            }
          }
          /* Old version
          for(auto nu : function_range(obs, dfbs_)) {
            for(auto sigma : function_range(obsptr, dfbsptr, 0, nu)) {
              //----------------------------------------//
              // Get the coefficients
              IntPair nu_sigma(nu, sigma);
              auto cpair = coefs_[nu_sigma];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              //----------------------------------------//
              // Add in the contribution to Kt(mu, nu)
              for(auto mu : function_range(ish)) {
                Kt_part(mu, nu) +=
                    B_mus[mu.bfoff_in_shell].row(sigma).segment(nu.atom_dfbfoff, nu.atom_dfnbf) * Ca;
                if(nu != sigma){
                  Kt_part(mu, sigma) +=
                      B_mus[mu.bfoff_in_shell].row(nu).segment(nu.atom_dfbfoff, nu.atom_dfnbf) * Ca;
                }
                if(nu.center != sigma.center) {
                  Kt_part(mu, nu) +=
                      B_mus[mu.bfoff_in_shell].row(sigma).segment(sigma.atom_dfbfoff, sigma.atom_dfnbf) * Cb;
                  if(nu != sigma){
                    Kt_part(mu, sigma) +=
                        B_mus[mu.bfoff_in_shell].row(nu).segment(sigma.atom_dfbfoff, sigma.atom_dfnbf) * Cb;
                  }
                }
              } // end loop over mu
              //----------------------------------------//
              // Add in the contribution to Kt(rho, nu)
              if(ish != jsh){
                for(auto rho : function_range(jsh)) {
                  Kt_part(rho, nu) +=
                      B_rhos[rho.bfoff_in_shell].row(sigma).segment(nu.atom_dfbfoff, nu.atom_dfnbf) * Ca;
                  if(nu != sigma){
                    Kt_part(rho, sigma) +=
                        B_rhos[rho.bfoff_in_shell].row(nu).segment(nu.atom_dfbfoff, nu.atom_dfnbf) * Ca;
                  }
                  if(nu.center != sigma.center) {
                    Kt_part(rho, nu) +=
                        B_rhos[rho.bfoff_in_shell].row(sigma).segment(sigma.atom_dfbfoff, sigma.atom_dfnbf) * Cb;
                    if(nu != sigma){
                      Kt_part(rho, sigma) +=
                          B_rhos[rho.bfoff_in_shell].row(nu).segment(sigma.atom_dfbfoff, sigma.atom_dfnbf) * Cb;
                    }
                  } // end if nu.center != sigma.center
                } // end loop over rho
              } // end if ish != jsh
              //----------------------------------------//
            } // end loop over sigma
          } // end loop over nu
          */
          if(ithr == 0) timer.exit("K contributions");
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
        } // end while get_shelsh)
        #endif


////////////////////////////////////////////////////////////////////////////////
// Two body part of K

                    //for(auto Y : function_range(Yblk)) {
                    //  Ct_mus[mu.bfoff_in_shell](rho, Y.bfoff_in_block) += 2.0 *
                    //      g2.row(Y.bfoff_in_block) * Cmu.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                    //}
                    //for(auto X : function_range(Xblk)) {
                    //  Ct_mus[mu.bfoff_in_shell](rho, Y.bfoff_in_block) += 2.0 *
                    //      g2(Y.bfoff_in_block, X.bfoff_in_block) * Cmu(X.bfoff_in_atom);
                    //}

                    //for(auto Y : function_range(Yblk)) {
                    //  Ct_mus[mu.bfoff_in_shell](rho, Y.bfoff_in_block) += 2.0 *
                    //      g2.row(Y.bfoff_in_block) * Crho.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                    //}
                    //for(auto X : function_range(Xblk)) {
                    //  Ct_mus[mu.bfoff_in_shell](rho, Y.bfoff_in_block) += 2.0 *
                    //      g2(Y.bfoff_in_block, X.bfoff_in_block) * Crho(X.bfoff_in_atom);
                    //}

        #if OLD_K
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          for(auto Yblk : shell_block_range(dfbs_, 0, 0, NoLastIndex, NoRestrictions)) {
            //----------------------------------------//
            std::vector<Eigen::MatrixXd> Ct_mus(ish.nbf);
            for(auto mu : function_range(ish)) {
              Ct_mus[mu.bfoff_in_shell].resize(jsh.nbf, Yblk.nbf);
              Ct_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(jsh.nbf, Yblk.nbf);
            }
            /*-----------------------------------------------------*/
            /* Form Ct_mu for each mu in ish                  {{{2 */ #if 2 // begin fold
            if(ithr == 0) timer.enter("form Ct_mu");
            //for(auto Xsh : iter_shells_on_center(dfbs_, ish.center)) {
            for(auto Xblk : iter_shell_blocks_on_center(dfbs_, ish.center)) {
              if(ithr == 0) timer.enter("compute ints");
              auto g2_part_ptr = ints_to_eigen(
                Yblk, Xblk,
                eris_2c_[ithr],
                coulomb_oper_type_
              );
              if(ithr == 0) timer.exit("compute ints");
              const auto& g2_part = *g2_part_ptr;
              for(auto mu : function_range(ish)){
                for(auto rho : function_range(jsh)){
                  //----------------------------------------//
                  // Get the coefficients
                  IntPair mu_rho(mu, rho);
                  auto cpair = coefs_[mu_rho];
                  VectorMap& Ca = *(cpair.first);
                  //----------------------------------------//
                  Ct_mus[mu.bfoff_in_shell].row(rho.bfoff_in_shell).transpose() +=
                      g2_part * Ca.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                  //----------------------------------------//
                } // end loop over functions rho in jsh
              } // end loop over functions mu in ish
            } // end loop of Xsh over shells on ish.center
            //----------------------------------------//
            if(ish.center != jsh.center){
              for(auto Xblk : iter_shell_blocks_on_center(dfbs_, jsh.center)) {
                auto g2_part_ptr = ints_to_eigen(
                  Yblk, Xblk,
                  eris_2c_[ithr],
                  coulomb_oper_type_
                );
                const auto& g2_part = *g2_part_ptr;
                for(auto mu : function_range(ish)){
                  for(auto rho : function_range(jsh)){
                    //----------------------------------------//
                    // Get the coefficients
                    IntPair mu_rho(mu, rho);
                    auto cpair = coefs_[mu_rho];
                    VectorMap& Cb = *(cpair.second);
                    //----------------------------------------//
                    Ct_mus[mu.bfoff_in_shell].row(rho.bfoff_in_shell).transpose() +=
                        g2_part * Cb.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                    //----------------------------------------//
                  } // end loop over functions rho in jsh
                } // end loop over functions mu in ish
              } // end loop of Xsh over shells on jsh.center
            } // end if ish.center != jsh.center
            //----------------------------------------//
            /********************************************************/ #endif //2}}}
            /*-----------------------------------------------------*/
            /* Form dt_mu, dt_rho for functions in ish, jsh   {{{2 */ #if 2 // begin fold
            if(ithr == 0) timer.change("form dt_mu, dt_rho");
            std::vector<Eigen::MatrixXd> dt_mus(ish.nbf);
            std::vector<Eigen::MatrixXd> dt_rhos((ish == jsh) ? 0 : (int)jsh.nbf);
            for(auto mu : function_range(ish)) {
              dt_mus[mu.bfoff_in_shell].resize(nbf, Yblk.nbf);
              dt_mus[mu.bfoff_in_shell] = 2.0 * D.middleCols(jsh.bfoff, jsh.nbf) * Ct_mus[mu.bfoff_in_shell];
              /*
              if(xml_debug) {
                Eigen::MatrixXd tmp(nbf, dfnbf);
                tmp = Eigen::MatrixXd::Zero(nbf, dfnbf);
                tmp.middleCols(Ysh.bfoff, Ysh.nbf) = dt_mus[mu.bfoff_in_shell];
                write_as_xml(
                    "new_dt_part", tmp,
                    std::map<std::string, int>{ {"mu", mu}, {"jsh", jsh} }
                );
              }
              */
            }
            if(ish != jsh){
              for(auto rho : function_range(jsh)) {
                dt_rhos[rho.bfoff_in_shell].resize(nbf, Yblk.nbf);
                dt_rhos[rho.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, Yblk.nbf);
                for(auto mu : function_range(ish)) {
                  dt_rhos[rho.bfoff_in_shell] += 2.0 * D.col(mu) * Ct_mus[mu.bfoff_in_shell].row(rho.bfoff_in_shell);
                }
                /*
                if(xml_debug) {
                  Eigen::MatrixXd tmp(nbf, dfnbf);
                  tmp = Eigen::MatrixXd::Zero(nbf, dfnbf);
                  tmp.middleCols(Ysh.bfoff, Ysh.nbf) = dt_rhos[rho.bfoff_in_shell];
                  write_as_xml(
                      "new_dt_part", tmp,
                      std::map<std::string, int>{ {"mu", rho}, {"jsh", ish} }
                  );
                }
                */
              }
            }
            /********************************************************/ #endif //2}}}
            /*-----------------------------------------------------*/
            /* Add contributions to Kt_part                   {{{2 */ #if 2 // begin fold
            if(ithr == 0) timer.change("K contributions");
            //----------------------------------------//
            for(auto Y : function_range(Yblk)) {
              Eigen::MatrixXd& C_Y = coefs_transpose_[Y];
              const int obs_atom_bfoff = obs->shell_to_function(obs->shell_on_center(Y.center, 0));
              const int obs_atom_nbf = obs->nbasis_on_center(Y.center);
              for(auto mu : function_range(ish)) {
                // dt_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
                // C_Y is (Y.{obs_}atom_nbf x nbf)
                // result should be (Y.{obs_}atom_nbf x 1)
                Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                    0.5 * C_Y * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block);
                // The sigma <-> nu term
                Kt_part.row(mu).transpose() -= 0.5
                    * C_Y.transpose()
                    * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
                // Add back in the nu.center == Y.center part
                Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() += 0.5
                    * C_Y.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                    * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
                //----------------------------------------//
                /* Failsafe version:
                //for(auto nu : iter_functions_on_center(obs, Y.center)) {
                //  //for(auto sigma : function_range(obs)) {
                //  //  Kt_part(mu, nu) -= 0.5 * C_Y(nu.bfoff_in_atom, sigma) * dt_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                //  //}
                //  Kt_part(mu, nu) -= 0.5 * C_Y.row(nu.bfoff_in_atom) * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell);
                //}
                //----------------------------------------//
                //for(auto nu : function_range(obs)) {
                //  if(nu.center != Y.center) {
                //    //for(auto sigma : iter_functions_on_center(obs, Y.center)) {
                //    //  Kt_part(mu, nu) -= 0.5 * C_Y(sigma.bfoff_in_atom, nu) * dt_mus[mu.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                //    //}
                //    Kt_part(mu, nu) -= 0.5 * C_Y.col(nu).transpose() * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_shell).segment(
                //      obs->shell_to_function(obs->shell_on_center(Y.center, 0)),
                //      obs->nbasis_on_center(Y.center)
                //    );
                //  }
                //}
                */
                //----------------------------------------//
              }
              if(ish != jsh){
                for(auto rho : function_range(jsh)) {
                  // dt_rhos[rho.bfoff_in_shell] is (nbf x Ysh.nbf)
                  // C_Y is (Y.{obs_}atom_nbf x nbf)
                  // result should be (Y.{obs_}atom_nbf x 1)
                  Kt_part.row(rho).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                      0.5 * C_Y * dt_rhos[rho.bfoff_in_shell].col(Y.bfoff_in_block);
                  // The sigma <-> nu term
                  Kt_part.row(rho).transpose() -= 0.5
                      * C_Y.transpose()
                      * dt_rhos[rho.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
                  // Add back in the nu.center == Y.center part
                  Kt_part.row(rho).segment(obs_atom_bfoff, obs_atom_nbf).transpose() += 0.5
                      * C_Y.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                      * dt_rhos[rho.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
                  //----------------------------------------//
                  /* Failsafe version:
                  //for(auto nu : function_range(obs)) {
                  //  for(auto sigma : function_range(obs)) {
                  //    if(nu.center == Y.center){
                  //      Kt_part(rho, nu) -= 0.5 * C_Y(nu.bfoff_in_atom, sigma) * dt_rhos[rho.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                  //    }
                  //    else if(sigma.center == Y.center){
                  //      Kt_part(rho, nu) -= 0.5 * C_Y(sigma.bfoff_in_atom, nu) * dt_rhos[rho.bfoff_in_shell](sigma, Y.bfoff_in_shell);
                  //    }
                  //  }
                  //}
                  */
                  //----------------------------------------//
                }
              }
            }
            //----------------------------------------//
            /* Old version
            for(auto ksh : shell_range(obs, dfbs_)) {
              for(auto lsh : shell_range(obs, dfbs_, 0, ksh)) {
              //for(auto lsh : iter_significant_partners(ksh)) {
                if(ksh.center != Ysh.center and lsh.center != Ysh.center) continue;
                //----------------------------------------//
                //const double pf = ksh == lsh ? 1.0 : 2.0;
                //----------------------------------------//
                //for(int nu = ksh.bfoff; nu <= ksh.last_function; ++nu){ // function_range(ksh)) {
                //  for(int sigma = lsh.bfoff; sigma <= lsh.last_function; ++sigma){ // function_range(ksh)) {
                for(auto nu : function_range(ksh)) {
                  for(auto sigma : function_range(lsh)) {
                    //----------------------------------------//
                    // Get the coefficients
                    IntPair nu_sigma(nu, sigma);
                    auto cpair = coefs_[nu_sigma];
                    VectorMap& C = (Ysh.center == ksh.center) ? *(cpair.first) : *(cpair.second);
                    //----------------------------------------//
                    //for(auto Y : function_range(Ysh)) {
                    //  auto& C_Y = coefs_transpose_[Y];
                    //  for(auto mu : function_range(ish)) {
                    //    Kt_part(mu, nu) -= 0.5
                    //        * dt_mus[mu.bfoff_in_shell].row(sigma)
                    //        * C_Y(nu.bfoff_in_atom, sigma);
                    //  }
                    //}
                    //----------------------------------------//
                    for(auto mu : function_range(ish)) {
                      Kt_part(mu, nu) -= 0.5
                          * dt_mus[mu.bfoff_in_shell].row(sigma)
                          * C.segment(Ysh.bfoff_in_atom, Ysh.nbf);
                      if(ksh != lsh){
                        Kt_part(mu, sigma) -= 0.5
                            * dt_mus[mu.bfoff_in_shell].row(nu)
                            * C.segment(Ysh.bfoff_in_atom, Ysh.nbf);
                      }
                    } // end loop over mu in ish
                    //----------------------------------------//
                    if(ish != jsh){
                      for(auto rho : function_range(jsh)) {
                        Kt_part(rho, nu) -= 0.5
                            * dt_rhos[rho.bfoff_in_shell].row(sigma)
                            * C.segment(Ysh.bfoff_in_atom, Ysh.nbf);
                        if(ksh != lsh){
                          Kt_part(rho, sigma) -= 0.5
                              * dt_rhos[rho.bfoff_in_shell].row(nu)
                              * C.segment(Ysh.bfoff_in_atom, Ysh.nbf);
                        }
                      } // end loop over rho in jsh
                    } // end if ish != jsh
                    //----------------------------------------//
                  } // end loop over sigma in lsh
                } // end loop over nu in ksh
                //----------------------------------------//
              } // end loop over lsh
            } // end loop over ksh
            */
            //----------------------------------------//
            if(ithr == 0) timer.exit();
            /********************************************************/ #endif //2}}}
            /*-----------------------------------------------------*/
          } // end loop over shells Ysh
        } // end while get_shell_pair(ish, jsh)
        #endif

////////////////////////////////////////////////////////////////////////////////
// other stuff

    /*
    shell_iter_arbitrary_wrapper<std::vector<int>>
    iter_significant_partners(
        const ShellData& ish
    )
    {
      return shell_iter_arbitrary_wrapper<std::vector<int>>(
          shell_to_sig_shells_[ish.index],
          ish.basis,
          ish.dfbasis
      );
    }
    */


// From J
            /*
            for(auto mu : function_range(ish)){
              for(auto nu : function_range(jsh)){
              //for(int jbf = 0, nu = jsh.bfoff; jbf < jsh.nbf; ++jbf, ++nu){
                //----------------------------------------//
                const int ijbf = mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell;
                //----------------------------------------//
                // compute dtilde contribution
                //dt.segment(Xblk.bfoff, Xblk.nbf) += perm_fact * D(mu, nu) * g_part->row(ijbf);
                //----------------------------------------//
                // add C_tilde * g_part to J
                if(nu <= mu){
                  // only constructing the lower triangle of the J matrix
                  jpart(mu, nu) += g_part->row(ijbf) * C_tilde.segment(Xblk.bfoff, Xblk.nbf);
                }
                //----------------------------------------//
              } // end loop over jbf
            } // end loop over ibf
            */
          //for(auto kshdf : shell_range(dfbs_)){
          //  std::shared_ptr<Eigen::MatrixXd> g_part = ints_to_eigen(
          //      ish, jsh, kshdf,
          //      eris_3c_[ithr],
          //      coulomb_oper_type_
          //  );
          //  for(int ibf = 0, mu = ish.bfoff; ibf < ish.nbf; ++ibf, ++mu){
          //    //for(auto jbf : function_range(jsh)){
          //    for(int jbf = 0, nu = jsh.bfoff; jbf < jsh.nbf; ++jbf, ++nu){
          //      //----------------------------------------//
          //      const int ijbf = ibf*jsh.nbf + jbf;
          //      //----------------------------------------//
          //      // compute dtilde contribution
          //      dt.segment(kshdf.bfoff, kshdf.nbf) += perm_fact * D(mu, nu) * g_part->row(ijbf);
          //      //----------------------------------------//
          //      // add C_tilde * g_part to J
          //      if(nu <= mu){
          //        // only constructing the lower triangle of the J matrix
          //        jpart(mu, nu) += g_part->row(ijbf) * C_tilde.segment(kshdf.bfoff, kshdf.nbf);
          //      }
          //      //----------------------------------------//
          //    } // end loop over jbf
          //  } // end loop over ibf
          //} // end loop over kshbf

#endif /* CADF_ATTIC_H_ */

