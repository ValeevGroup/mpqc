//
// assignments.cc
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

#include "cadfclhf.h"
#include "assignments.h"

using namespace sc;
using namespace sc::cadf;


AssignmentGrid::AssignmentGrid(
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis,
    int n_node
) : basis_(basis),
    dfbasis_(dfbasis)
{
  for(auto&& iblk : shell_block_range(basis, dfbasis, 0, NoLastIndex, SameCenter, NoMaximumBlockSize)) {
    assert(iblk.nbf == iblk.atom_nbf);
    atoms_.emplace_back(iblk);
  }
  for(auto&& ish : shell_range(basis, dfbasis)) {
    obs_shells_.emplace_back(ish);
  }
  const int natoms = atoms_.size();

  // We want the number of df rows to always be less than or equal to the number of atoms
  const int nrows_dfbs = floor(sqrt((double)n_node));
  const int nrows_obs = n_node / nrows_dfbs;

  for(int irow = 0; irow < nrows_obs; ++irow) {
    auto handle = obs_rows_.emplace(irow, false);
    (*handle).pq_handle = handle;
  }

  for(int irow = 0; irow < nrows_dfbs; ++irow) {
    auto handle = dfbs_rows_.emplace(irow, true);
    (*handle).pq_handle = handle;
  }

  for(auto& atom : atoms_) {
    // Assign atom to row with smallest workload
    const auto& dfbs_row = dfbs_rows_.top();
    (*dfbs_row.pq_handle).assign_item(&atom);
    dfbs_rows_.update(dfbs_row.pq_handle);
  }

  for(auto& shell : obs_shells_) {
    // And the same for the obs rows
    const auto& obs_row = obs_rows_.top();
    (*obs_row.pq_handle).assign_item(&shell);
    obs_rows_.update(obs_row.pq_handle);
  }

  uint bin_id = 0;
  for(auto&& dfbsrow : dfbs_rows_) {
    for(auto&& obsrow : obs_rows_) {
      auto handle = bins_.emplace(
          bin_id++, this
      );
      (*handle).pq_handle = handle;
      (*handle).register_in_row(obsrow, false);
      (*handle).register_in_row(dfbsrow, true);
    }
  }

  for(auto&& row : dfbs_rows_) {
    (*row.pq_handle).make_assignments();
  }
  for(auto&& row : obs_rows_) {
    (*row.pq_handle).make_assignments();
  }

  for(int inode = 0; inode < n_node; ++inode) {
    const AssignmentBin& most_work_bin = bins_.top();
    nodes_.push_back((*most_work_bin.pq_handle).add_node(inode));
    bins_.update(most_work_bin.pq_handle);
  }

  for(auto&& bin : bins_) {
    (*bin.pq_handle).make_assignments();
  }

}

void
AssignmentGrid::print_detail(std::ostream& o) const
{
  o << indent << "Total number of bins: " << bins_.size()
    << " (" << obs_rows_.size() << " obs rows by " << dfbs_rows_.size() << " dfbs rows)"
    << std::endl;
  o << indent << "Node-wise coefficient memory requirements:" << std::endl;
  const int nnode = nodes_.size();

  o << incindent << indent;
  constexpr int n_per_row = 6;
  for(auto node : nodes_) {
    o << std::setw(12) << std::right << data_size_to_string((node->bin->obs_ncoefs + node->bin->dfbs_ncoefs)*sizeof(double));
    if(node->id + 1 != nnode && (node->id + 1) % n_per_row == 0) {
      o << std::endl << indent;
    }
  }
  o << decindent << std::endl;

  o << indent << "Shell pairs to compute per node:" << std::endl;
  o << incindent << indent;
  for(auto node : nodes_) {
    o << std::setw(12) << std::right << node->shell_pair_count;
    if(node->id + 1 != nnode && (node->id + 1) % n_per_row == 0) {
      o << std::endl << indent;
    }
  }

  o << decindent << std::endl;
  o << indent << "Basis function pairs to compute per node:" << std::endl;
  o << incindent << indent;
  for(auto node : nodes_) {
    o << std::setw(12) << std::right << node->basis_pair_count;
    if(node->id + 1 != nnode && (node->id + 1) % n_per_row == 0) {
      o << std::endl << indent;
    }
  }
  o << decindent << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void
AssignmentBinRow::make_assignments()
{
  for(auto item_ptr : assigned_items) {
    auto bin_handle =  bins.top()->row_handles[is_df_row];
    (*bin_handle)->compute_coef_for_item(item_ptr, is_df_row);
    bins.update(bin_handle);
  }
}

////////////////////////////////////////////////////////////////////////////////

AssignmentBin::AssignmentBin(
    uint id,
    AssignmentGrid* grid
) : id(id), grid(grid), obs_row_id(-1), dfbs_row_id(-1)
{

}

bool
AssignmentBin::operator<(const AssignmentBin& other) const
{
  // We want to assign nodes to bins first that don't have any nodes and then to the bin with
  // the greatest workload per node; the priority queue, which is a max heap, should be ordered in this manor
  if(n_node() == 0) {
    if(other.n_node() == 0) return estimated_workload < other.estimated_workload;
    else return false; // in other words, this should have higher priority than other
  }
  else if(other.n_node() == 0) return true; // other should have higher priority than this
  else {
    return double(estimated_workload)/double(std::max(n_node(), (size_t)1))
        > double(other.estimated_workload)/double(std::max(other.n_node(), (size_t)1));
  }
}

void
AssignmentBin::register_in_row(const AssignmentBinRow& row, bool is_df) {
  (*row.pq_handle).add_bin(this);
  if(is_df) {
    dfbs_row_id = row.id;
    for(auto& dfbs_atom : row.assigned_items) {
      assign_dfbs_atom(dfbs_atom);
    }
  }
  else {
    obs_row_id = row.id;
    for(auto& obs_shell : row.assigned_items) {
      assign_obs_shell(obs_shell);
    }
  }
}

cadf::Node*
AssignmentBin::add_node(int index)
{
  nodes_list.emplace_back();
  nodes_list.back().node_index = index;
  nodes_list.back().bin = this;
  auto handle = nodes.push(&nodes_list.back());
  (*handle)->pq_handle = handle;
  return *handle;
}

inline void
AssignmentBin::make_assignments()
{
  for(auto dfbs_atom : compute_coef_items[true]) {
    auto handle = nodes.top()->pq_handle;
    (*handle)->assign_coef_item(dfbs_atom, true);
    nodes.update(handle);
  }
  for(auto obs_shell : compute_coef_items[false]) {
    auto handle = nodes.top()->pq_handle;
    (*handle)->assign_coef_item(obs_shell, false);
    nodes.update(handle);
  }

  for(auto node : nodes) {
    (*node->pq_handle)->estimated_workload = 0;
    nodes.update(node->pq_handle);
  }

  assert(n_node() > 0);
  for(auto dfbs_atom : assigned_dfbs_atoms) {
    for(auto obs_shell : assigned_obs_shells) {
      auto Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, grid->dfbasis(), grid->basis());
      auto ish = ShellData(obs_shell->index, grid->basis(), grid->dfbasis());
      auto& top_node = nodes.top();
      (*top_node->pq_handle)->assign_pair(ish, Xblk);
      nodes.update(top_node->pq_handle);
    }
  }

  // Update the offsets
  for(auto& off_pair : dfbs_coef_offsets) {
    off_pair.second += obs_ncoefs;
  }

}
