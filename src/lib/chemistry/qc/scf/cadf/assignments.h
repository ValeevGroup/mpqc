//
// assignments.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Apr 15, 2014
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

#ifndef _chemistry_qc_scf_assignments_h
#define _chemistry_qc_scf_assignments_h

#include <array>
#include <set>
#include <utility>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
//#include <boost/container/set.hpp>


#include <boost/heap/fibonacci_heap.hpp>

#include <util/misc/formio.h>

#include "iters.h"


namespace sc { namespace cadf {

namespace detail {
  template<typename T, template<typename...> class compare=std::less>
  struct deref_compare;

  template<typename T, template <typename...> class compare>
  struct deref_compare<T*, compare>
  {
      bool operator()(T* const& a, T* const& b) const {
        return compare<T>()(*a, *b);
      }
    private:
      //compare<T> cmp_;
  };

  template<typename T, template <typename...> class compare>
  struct deref_compare<boost::shared_ptr<T>, compare>
  {
      bool operator()(boost::shared_ptr<T> const& a, boost::shared_ptr<T> const& b) const {
        return compare<T>()(*a, *b);
      }
    private:
      //compare<T> cmp_;
  };


  template<typename T>
  struct more_work {
      bool operator()(const T& a, const T& b) const {
        return a.coef_workload > b.coef_workload;
      }
  };

  template<typename T>
  struct index_less {
      bool operator()(const T& a, const T& b) const {
        return a.index < b.index;
      }
  };

}

template<typename T, typename... Args> using priority_queue =
    boost::heap::fibonacci_heap<T, boost::heap::stable<true>, Args...>;
template<typename T, template<typename...> class compare=std::less> using ptr_priority_queue =
    boost::heap::fibonacci_heap<T, boost::heap::stable<true>,
      boost::heap::compare<
        detail::deref_compare<T, compare>
      >
    >;
template<
  template<typename...> class container,
  typename T,
  template<typename...> class compare=std::less
> using ordered_ptr_container =
    container<
      T, detail::deref_compare<T, compare>
    >;
template<typename T, template<typename...> class compare=std::less> using ptr_set =
    ordered_ptr_container<std::set, T, compare>;
typedef uint64_t uli;
typedef unsigned int uint;

class Node;

class AssignableItem {
  public:
    int index;
    uli coefs_size;
    AssignableItem(int index) : index(index), coefs_size(0) { }
    virtual ~AssignableItem() { }
    /// Cost estimate for atom's portion of the coefficient tensor
    virtual uli cost_estimate(bool df) const = 0;
};

class AssignableAtom : public AssignableItem {
  public:
    int nshell;
    int nbf;
    int dfnshell;
    int dfnbf;

    // Note: coefs_size is the size of coefs
    //   C_{\nu_c \sigma}^{X_c} for atom c (in number of doubles)

    AssignableAtom(const ShellBlockData<>& iblk)
      : AssignableItem(iblk.center)
    {
      assert(iblk.nbf == iblk.atom_nbf);
      nshell = iblk.nshell;
      nbf = iblk.nbf;
      dfnshell = iblk.atom_dfnsh;
      dfnbf = iblk.atom_dfnbf;
      coefs_size = iblk.nbf * iblk.basis->nbasis() * iblk.atom_dfnbf;
    }

    uli cost_estimate(bool df) const {
      return coefs_size;
    }
};

class AssignableShell : public AssignableItem {
  public:

    int nbf;

    // Note: coefs_size is the size of coefs
    //   C_{\mu_a \sigma}^{X_a} for shell \mu_a (in number of doubles)

    AssignableShell(const ShellData& ish)
      : AssignableItem(ish.index)
    {
      nbf = ish.nbf;
      coefs_size = ish.nbf * ish.basis->nbasis() * ish.atom_dfnbf;
    }

    uli cost_estimate(bool df) const {
      return coefs_size;
    }
};

class AssignableShellPair {

    uli cost_estimate_;

  public:
    int ish;
    int Xatom;
    const boost::shared_ptr<Node> node;

    AssignableShellPair(
        const ShellData& ish,
        const ShellBlockData<>& Xblk,
        const boost::shared_ptr<Node> node
    ) : ish(ish), Xatom(Xblk.center), node(node)
    {
      // TODO more efficient cost estimation options
      cost_estimate_ = ish.nbf * Xblk.nbf;
    }

    uli cost_estimate() const {
      return cost_estimate_;
    }
};


class AssignmentBin;
namespace assignments { struct AtomCluster; };

class Node : public boost::enable_shared_from_this<Node> {
  public:

    union { int node_index; int id; };
    std::vector<AssignableShellPair> pairs;
    std::set<uli> obs_shells_to_do;
    std::set<uli> dfbs_atoms_to_do;
    std::array<std::vector<boost::shared_ptr<AssignableItem>>, 2> compute_coef_items;
    boost::shared_ptr<AssignmentBin> bin;
    boost::shared_ptr<assignments::AtomCluster> cluster;
    typename ptr_priority_queue<boost::shared_ptr<Node>>::handle_type pq_handle;
    uli estimated_workload = 0;
    uli shell_pair_count = 0;
    uli basis_pair_count = 0;

    bool is_me = false;

    uli assign_pair(
        const ShellData& ish,
        const ShellBlockData<>& Xblk
    )
    {
      obs_shells_to_do.insert(ish);
      dfbs_atoms_to_do.insert(Xblk.center);
      AssignableShellPair pair(ish, Xblk, shared_from_this());
      if(is_me) {
        pairs.emplace_back(pair);
      }
      // compute the cost
      const uli cost = pair.cost_estimate();
      estimated_workload += cost;
      shell_pair_count += Xblk.nshell;
      basis_pair_count += Xblk.nbf * ish.nbf;
      return cost;
    }

    bool should_do_obs_shell(uli shell_index) const {
      return obs_shells_to_do.find(shell_index) != obs_shells_to_do.end();
    }

    bool should_do_dfbs_atom(uli atom_index) const {
      return dfbs_atoms_to_do.find(atom_index) != dfbs_atoms_to_do.end();
    }

    void assign_coef_item(boost::shared_ptr<AssignableItem> const& item, bool is_df) {
      compute_coef_items[is_df].push_back(item);
      estimated_workload += item->cost_estimate(is_df);
    }

    uli dfnsh() const;

    const std::set<uint>& assigned_dfbs_shells() const;

    const ptr_set<boost::shared_ptr<AssignableAtom>, detail::index_less>&
    assigned_dfbs_atoms() const;



    bool operator <(const Node& other) const {
      return estimated_workload > other.estimated_workload;
    }

    const std::unordered_map<uint, uli>&
    dfbs_coef_offsets() const;

};

class AssignmentGrid;
class AssignmentBinRow;

class AssignmentBin : public boost::enable_shared_from_this<AssignmentBin> {


  public:
    ptr_priority_queue<boost::shared_ptr<Node>> nodes;
    std::vector<boost::shared_ptr<Node>> nodes_list;
    ptr_set<boost::shared_ptr<AssignableAtom>, detail::index_less> assigned_dfbs_atoms;
    std::set<uint> assigned_dfbs_shells;
    ptr_set<boost::shared_ptr<AssignableItem>, detail::index_less> assigned_obs_shells;
    std::array<std::vector<boost::shared_ptr<AssignableItem>>, 2> compute_coef_items;
    uli estimated_workload = 0;
    uli coef_workload = 0;
    uint id;
    uint obs_row_id;
    uint dfbs_row_id;

    // Size of coefs C_{\mu_a \rho}^{X_a} for a in assigned obs_atoms
    uli obs_ncoefs = 0;
    // Size of coefs C_{\nu_c \sigma}^{X_c} for c in assigned dfbs_atoms
    uli dfbs_ncoefs = 0;
    // offsets of the start of C_{\mu_a \rho}^{X_a} in the obs coefficient memory (in number of doubles) for given a
    std::unordered_map<uint, uli> obs_coef_offsets;
    // offsets of the start of C_{\nu_c \rho}^{X_c} in the dfbs coefficient memory (in number of doubles) for given c
    std::unordered_map<uint, uli> dfbs_coef_offsets;

    AssignmentGrid* grid;

    AssignmentBin(
        uint id, AssignmentGrid* grid
    );

    boost::shared_ptr<cadf::Node> add_node(int index);

    inline void register_in_row(const AssignmentBinRow& row, bool is_df);

    inline void assign_dfbs_atom(const boost::shared_ptr<AssignableItem>& dfbs_atom);

    inline void assign_obs_shell(const boost::shared_ptr<AssignableItem>& obs_shell);

    inline void make_assignments();

    void compute_coef_for_item(const boost::shared_ptr<AssignableItem>& item, bool is_df);

    size_t n_node() const { return nodes.size(); }

    uint dfnsh() const { return assigned_dfbs_shells.size(); }

    bool operator<(const AssignmentBin& other) const;

    typename ptr_priority_queue<boost::shared_ptr<AssignmentBin>>::handle_type pq_handle;
    std::array<typename ptr_priority_queue<boost::shared_ptr<AssignmentBin>, detail::more_work>::handle_type, 2> row_handles;

  private:
    static constexpr bool debug_ = false;
};

class AssignmentBinRow {
  public:
    bool is_df_row;
    uli estimated_workload = 0;
    std::vector<boost::shared_ptr<AssignableItem>> assigned_items;
    ptr_priority_queue<boost::shared_ptr<AssignmentBin>, detail::more_work> bins;
    uint id;

    AssignmentBinRow(uint id, bool is_df)
      : id(id), is_df_row(is_df)
    { }

    void assign_item(const boost::shared_ptr<AssignableItem>& item) {
      assigned_items.push_back(item);
      estimated_workload += item->cost_estimate(is_df_row);
    }

    void add_bin(boost::shared_ptr<AssignmentBin> const& bin) {
      auto handle = bins.push(bin);
      (*handle)->row_handles[is_df_row] = handle;
    }

    void make_assignments();

    bool operator <(const AssignmentBinRow& other) const {
      // We want to assign work to the row with the lowest estimated workload
      return estimated_workload > other.estimated_workload;
    }

    typename priority_queue<AssignmentBinRow>::handle_type pq_handle;
};

class AssignmentGrid {

    GaussianBasisSet* basis_;
    GaussianBasisSet* dfbasis_;

    // TODO Fix this.  We can't store pointers to these things in other places, since they may be moved
    std::vector<boost::shared_ptr<AssignableAtom>> atoms_;
    std::vector<boost::shared_ptr<AssignableShell>> obs_shells_;

    ptr_priority_queue<boost::shared_ptr<AssignmentBin>> bins_;
    priority_queue<AssignmentBinRow> obs_rows_;
    priority_queue<AssignmentBinRow> dfbs_rows_;
    std::vector<boost::shared_ptr<Node>> nodes_;

  public:

    int nrows_dfbs;
    int nrows_obs;
    int nbin;
    bool bins_have_multiple_nodes = false;


    AssignmentGrid(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int n_node, int me
    );

    const AssignmentBin& my_bin(int me) const {
      return *(nodes_[me]->bin);
    }

    const Node& my_assignments(int me) const {
      return *(nodes_[me]);
    }

    const boost::shared_ptr<Node> my_assignments_ptr(int me) const {
      return nodes_[me];
    }

    GaussianBasisSet* basis() {
      return basis_;
    }

    GaussianBasisSet* dfbasis() {
      return dfbasis_;
    }

    void print_detail(std::ostream& o=ExEnv::out0(), bool full_memory=false) const;

};


namespace assignments {

struct Assignments;

struct AtomCluster : boost::enable_shared_from_this<AtomCluster> {
    ptr_set<boost::shared_ptr<AssignableAtom>, detail::index_less> atoms;
    std::vector<boost::shared_ptr<Node>> nodes;
    Assignments* parent;
    uli coefs_size = 0;
    int index;
    uli dfnsh = 0;
    std::set<uint> assigned_dfbs_shells;
    std::unordered_map<uint, uint> dfbs_shell_map;
    std::unordered_map<uint, uli> coef_offsets;

    void assign_atom(const boost::shared_ptr<AssignableAtom>& atom);

    boost::shared_ptr<Node> use_node(int inode, bool is_me) {
      nodes.emplace_back(boost::make_shared<Node>());
      auto& node = *nodes.back();
      node.id = inode;
      node.is_me = is_me;
      node.cluster = shared_from_this();
      return nodes.back();
    }

    uli workload_per_node() { return nodes.size() > 0 ? coefs_size / nodes.size() : std::numeric_limits<uli>::max(); }

    void make_assignments();

    bool has_local_coefs_for_atom(uint atom_index) {
      return coef_offsets.find(atom_index) != coef_offsets.end();
    }

};


struct Assignments {
    std::vector<boost::shared_ptr<AtomCluster>> clusters;
    std::vector<boost::shared_ptr<Node>> nodes;

    GaussianBasisSet* basis;
    GaussianBasisSet* dfbasis;

    Assignments(
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis,
        int min_atoms_per_cluster,
        int n_node, int me
    );

    void print_detail(std::ostream& o) const;


};

}

}} // end namespaces cadf and sc



#endif /* _chemistry_qc_scf_assignments_h */
