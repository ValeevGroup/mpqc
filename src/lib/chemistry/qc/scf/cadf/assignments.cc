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

using namespace std;


AssignmentGrid::AssignmentGrid(
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis,
    int n_node, int me
) : basis_(basis),
    dfbasis_(dfbasis)
{

  ExEnv::out0() << incindent;

  ExEnv::out0() << indent << "Filling list of dfbs atoms to assign" << endl;
  for(auto&& iblk : shell_block_range(basis, dfbasis, 0, NoLastIndex, SameCenter, NoMaximumBlockSize)) {
    assert(iblk.nbf == iblk.atom_nbf);
    atoms_.push_back(boost::make_shared<AssignableAtom>(iblk));
  }

  ExEnv::out0() << indent << "Filling list of obs shells to assign" << endl;
  for(auto&& ish : shell_range(basis, dfbasis)) {
    obs_shells_.push_back(boost::make_shared<AssignableShell>(ish));
  }
  const int natoms = atoms_.size();

  // We want the number of df rows to always be less than or equal to the number of atoms
  nrows_dfbs = floor(sqrt((double)n_node));
  nrows_obs = n_node / nrows_dfbs;
  nbin = nrows_dfbs * nrows_obs;

  ExEnv::out0() << indent << "Creating " << nrows_obs << " obs AssignmentBinRow objects" << endl;
  for(int irow = 0; irow < nrows_obs; ++irow) {
    // id is used for making MPI groups, so it must be unique.  Thus, offset it by number of bins
    //auto handle = obs_rows_.emplace(nbin + irow, false);
    auto handle = obs_rows_.push(AssignmentBinRow(nbin + irow, false));
    (*handle).pq_handle = handle;
  }

  ExEnv::out0() << indent << "Creating " << nrows_dfbs << " dfbs AssignmentBinRow objects" << endl;
  for(int irow = 0; irow < nrows_dfbs; ++irow) {
    // id is used for making MPI groups, so it must be unique.  Thus, offset it by number of bins + number of obs rows
    //auto handle = dfbs_rows_.emplace(nbin + nrows_obs + irow, true);
    auto handle = dfbs_rows_.push(AssignmentBinRow(nbin + nrows_obs + irow, true));
    (*handle).pq_handle = handle;
  }

  ExEnv::out0() << indent << "Assigning atoms to dfbs rows" << endl;
  for(auto& atom : atoms_) {
    // Assign atom to row with smallest workload
    const auto& dfbs_row = dfbs_rows_.top();
    (*dfbs_row.pq_handle).assign_item(atom);
    dfbs_rows_.update(dfbs_row.pq_handle);
  }

  ExEnv::out0() << indent << "Assigning shells to obs rows" << endl;
  for(auto& shell : obs_shells_) {
    // And the same for the obs rows
    const auto& obs_row = obs_rows_.top();
    (*obs_row.pq_handle).assign_item(shell);
    obs_rows_.update(obs_row.pq_handle);
  }

  ExEnv::out0() << indent << "Creating AssignmentBin objects" << endl;
  uint bin_id = 0;
  for(auto&& dfbsrow : dfbs_rows_) {
    for(auto&& obsrow : obs_rows_) {
      auto handle = bins_.push(boost::make_shared<AssignmentBin>(
          bin_id++, this
      ));
      (*handle)->pq_handle = handle;
      (*handle)->register_in_row(obsrow, false);
      (*handle)->register_in_row(dfbsrow, true);
    }
  }

  ExEnv::out0() << indent << "dfbs rows making assignments to bins" << endl;
  for(auto&& row : dfbs_rows_) {
    (*row.pq_handle).make_assignments();
  }
  ExEnv::out0() << indent << "obs rows making assignments to bins" << endl;
  for(auto&& row : obs_rows_) {
    (*row.pq_handle).make_assignments();
  }

  ExEnv::out0() << indent << "Assigning nodes to bins based on workload." << endl;
  for(int inode = 0; inode < n_node; ++inode) {
    const boost::shared_ptr<AssignmentBin>& most_work_bin = bins_.top();
    auto new_node = (*most_work_bin->pq_handle)->add_node(inode);
    if(inode == me) new_node->is_me = true;
    nodes_.push_back(new_node);
    bins_.update(most_work_bin->pq_handle);

  }

  ExEnv::out0() << indent << "AssignmentBins making assignments to available nodes" << endl;
  for(auto&& bin : bins_) {
    (*(bin->pq_handle))->make_assignments();
  }

  if(bins_have_multiple_nodes) {
    // TODO Fix this
    throw FeatureNotImplemented("More than one node per bin; use a perfect square number of nodes for now", __FILE__, __LINE__);
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Done making static distribution for exchange" << endl;

}

void
AssignmentGrid::print_detail(std::ostream& o, bool full_memory) const
{
  o << indent << "Total number of bins: " << bins_.size()
    << " (" << obs_rows_.size() << " obs rows by " << dfbs_rows_.size() << " dfbs rows)"
    << std::endl;
  o << indent << "Node-wise coefficient memory requirements:" << std::endl;
  const int nnode = nodes_.size();

  o << incindent << indent;
  constexpr int n_per_row = 6;
  for(auto node : nodes_) {
    if(full_memory) {
      o << std::setw(12) << std::right << "(all coefs)";
      if(node->id + 1 != nnode && (node->id + 1) % n_per_row == 0) {
        o << std::endl << indent;
      }

    }
    else {
      o << std::setw(12) << std::right << data_size_to_string((node->bin->obs_ncoefs + node->bin->dfbs_ncoefs)*sizeof(double));
      if(node->id + 1 != nnode && (node->id + 1) % n_per_row == 0) {
        o << std::endl << indent;
      }
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

uli
Node::dfnsh() const
{
  if(bin) {
    assert(!cluster);
    return bin->dfnsh();
  }
  else {
    assert(cluster);
    return cluster->dfnsh;
  }
}

const std::set<uint>&
Node::assigned_dfbs_shells() const
{
  if(bin) {
    assert(!cluster);
    return bin->assigned_dfbs_shells;
  }
  else {
    assert(cluster);
    return cluster->assigned_dfbs_shells;
  }
}

const ptr_set<boost::shared_ptr<AssignableAtom>, sc::cadf::detail::index_less>&
Node::assigned_dfbs_atoms() const
{
  if(bin) {
    assert(!cluster);
    return bin->assigned_dfbs_atoms;
  }
  else {
    assert(cluster);
    return cluster->atoms;
  }
}

const std::unordered_map<uint, uli>&
Node::dfbs_coef_offsets() const
{
  if(bin) {
    assert(!cluster);
    return bin->dfbs_coef_offsets;
  }
  else {
    assert(cluster);
    return cluster->coef_offsets;
  }
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
{ }

bool
AssignmentBin::operator<(const AssignmentBin& other) const
{
  // We want to assign nodes to bins first that don't have any nodes and then to the bin with
  // the greatest workload per node; the priority queue, which is a max heap, should be ordered in this manor
  if(n_node() == 0) {
    if(other.n_node() == 0) return estimated_workload > other.estimated_workload;
    else return false; // in other words, this should have higher priority than other
  }
  else if(other.n_node() == 0) return true; // other should have higher priority than this
  else {
    // TODO CHANGING THIS CHANGES THE ANSWER??!??!??
    return double(estimated_workload)/double(std::max(n_node(), (size_t)1))
        < double(other.estimated_workload)/double(std::max(other.n_node(), (size_t)1));
  }
}

inline void
AssignmentBin::register_in_row(const AssignmentBinRow& row, bool is_df) {
  (*row.pq_handle).add_bin(shared_from_this());
  if(is_df) {
    dfbs_row_id = row.id;
    for(auto dfbs_atom : row.assigned_items) {
      assign_dfbs_atom(dfbs_atom);
    }
  }
  else {
    obs_row_id = row.id;
    for(auto obs_shell : row.assigned_items) {
      assign_obs_shell(obs_shell);
    }
  }
}

boost::shared_ptr<cadf::Node>
AssignmentBin::add_node(int index)
{
  nodes_list.push_back(boost::make_shared<Node>());
  nodes_list.back()->node_index = index;
  nodes_list.back()->bin = shared_from_this();
  if(nodes_list.size() > 1) grid->bins_have_multiple_nodes = true;
  auto handle = nodes.push(nodes_list.back());
  (*handle)->pq_handle = handle;
  return *handle;
}

inline void
AssignmentBin::make_assignments()
{
  assert(n_node() > 0);

  if(debug_) {
    ExEnv::out0() << incindent;
    ExEnv::out0() << indent << "Making assignments in bin " << id << endl;
    ExEnv::out0() << incindent;
    ExEnv::out0() << indent << "Bin " << id << " has " << n_node() << " node" << (n_node() == 1 ? "" : "s")
                  << " to make assignments to." << endl;
  }

  if(debug_) ExEnv::out0() << indent << "Assigning coefficient computation for " << compute_coef_items[true].size() << " dfbs atoms" << endl;
  for(auto dfbs_atom : compute_coef_items[true]) {
    auto handle = nodes.top()->pq_handle;
    (*handle)->assign_coef_item(dfbs_atom, true);
    nodes.update(handle);
  }

  if(debug_) ExEnv::out0() << indent << "Assigning coefficient computation for " << compute_coef_items[false].size() << " obs shells" << endl;
  for(auto obs_shell : compute_coef_items[false]) {
    auto handle = nodes.top()->pq_handle;
    (*handle)->assign_coef_item(obs_shell, false);
    nodes.update(handle);
  }

  if(debug_) ExEnv::out0() << indent << "Resetting estimated workloads to make exchange matrix work assignments" << endl;
  for(auto node : nodes) {
    (*node->pq_handle)->estimated_workload = 0;
    nodes.update(node->pq_handle);
  }

  if(debug_) {
    ExEnv::out0() << indent << "Assigning pairs:" << endl;
    ExEnv::out0() << incindent;
    ExEnv::out0() << indent << assigned_dfbs_atoms.size() << " dfbs atoms, ";
    ExEnv::out0() << assigned_obs_shells.size() << " obs shells, ";
    ExEnv::out0() << assigned_dfbs_atoms.size() * assigned_obs_shells.size() << " total pairs" << endl;
    ExEnv::out0() << indent << "coef workload: " << coef_workload << ", exchange workload: " << estimated_workload << endl;
    ExEnv::out0() << indent << "estimated workload per node: " << double(estimated_workload) / double(n_node()) << endl;
    ExEnv::out0() << decindent;
  }
  for(auto dfbs_atom : assigned_dfbs_atoms) {
    for(auto obs_shell : assigned_obs_shells) {
      auto Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, grid->dfbasis(), grid->basis());
      auto ish = ShellData(obs_shell->index, grid->basis(), grid->dfbasis());
      auto& top_node = nodes.top();
      (*(top_node->pq_handle))->assign_pair(ish, Xblk);
      nodes.update(top_node->pq_handle);
    }
  }

  if(debug_) {
    ExEnv::out0() << indent << "Done making assignments for bin " << id << endl;
    ExEnv::out0() << decindent;
    ExEnv::out0() << decindent;
  }

}

inline void
AssignmentBin::assign_dfbs_atom(const boost::shared_ptr<AssignableItem>& dfbs_atom)
{
  assigned_dfbs_atoms.insert(boost::static_pointer_cast<AssignableAtom>(dfbs_atom));
  for(auto&& Xsh : iter_shells_on_center(grid->dfbasis(), dfbs_atom->index, grid->basis())) {
    assigned_dfbs_shells.insert(Xsh.index);
  }
  dfbs_coef_offsets[dfbs_atom->index] = dfbs_ncoefs;
  dfbs_ncoefs += dfbs_atom->coefs_size;
  estimated_workload += dfbs_ncoefs;
}

inline void
AssignmentBin::assign_obs_shell(const boost::shared_ptr<AssignableItem>& obs_shell)
{
  // MEMORY LEAK?
  //assigned_obs_shells.insert(boost::static_pointer_cast<AssignableShell>(obs_shell));
  assert(assigned_obs_shells.find(obs_shell) == assigned_obs_shells.end());
  assigned_obs_shells.insert(obs_shell);
  obs_coef_offsets[obs_shell->index] = obs_ncoefs;
  obs_ncoefs += obs_shell->coefs_size;
  estimated_workload += obs_ncoefs;
}

inline void
AssignmentBin::compute_coef_for_item(
    const boost::shared_ptr<AssignableItem>& item,
    bool is_df
)
{
  compute_coef_items[is_df].push_back(item);
  coef_workload += item->cost_estimate(is_df);
}

////////////////////////////////////////////////////////////////////////////////

using namespace sc::cadf::assignments;

Assignments::Assignments(
    GaussianBasisSet* basis,
    GaussianBasisSet* dfbasis,
    int min_atoms_per_cluster,
    int n_node, int me
) : basis(basis), dfbasis(dfbasis)
{
  const int natom = basis->ncenter();
  int ncluster = natom / min_atoms_per_cluster;
  ncluster = std::max(std::min(n_node, ncluster), 1);

  // Initialize clusters
  for(int i = 0; i < ncluster; ++i) {
    clusters.emplace_back(boost::make_shared<AtomCluster>());
    clusters.back()->index = i;
    clusters.back()->parent = this;
  }

  // Assign atoms to clusters
  // TODO better clustering algorithm that takes into account position as well as number of basis functions
  for(auto&& iblk : shell_block_range(basis, dfbasis, 0, NoLastIndex, SameCenter, NoMaximumBlockSize)) {
    assert(iblk.nbf == iblk.atom_nbf);
    clusters[0]->assign_atom(boost::make_shared<AssignableAtom>(iblk));
    int spot = 0;
    // Emulate a priority queue
    while(spot < clusters.size() - 1 and clusters[spot]->coefs_size > clusters[spot+1]->coefs_size) {
      std::swap(clusters[spot], clusters[spot+1]);
      ++spot;
    }
  }

  // Assign nodes to clusters based on workload
  for(int inode = 0; inode < n_node; ++inode) {
    int spot = clusters.size() - 1;
    nodes.push_back(clusters[spot]->use_node(inode, inode == me));
    while(spot > 0 and clusters[spot]->workload_per_node() < clusters[spot-1]->workload_per_node()) {
      std::swap(clusters[spot], clusters[spot-1]);
      --spot;
    }
  }

  // Let the clusters make assignements to it's nodes
  for(auto&& cluster : clusters) {
    cluster->make_assignments();
  }

}

void
AtomCluster::make_assignments()
{
  // Loop over shells
  for(auto&& ish : shell_range(parent->basis, parent->dfbasis)) {
    for(auto&& atom_ptr : atoms) {
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(atom_ptr->index, parent->basis, parent->dfbasis);
      nodes[0]->assign_pair(ish, Xblk);
    }
    int spot = 0;
    while(spot < nodes.size() - 1 and nodes[spot]->estimated_workload > nodes[spot+1]->estimated_workload) {
      std::swap(nodes[spot], nodes[spot+1]);
      ++spot;
    }
  }

  // Reset the workloads
  for(auto&& node : nodes) {
    node->estimated_workload = 0;
  }

  // Assign coefficients to compute
  for(auto atom_ptr : atoms) {
    ShellBlockData<> Xblk = ShellBlockData<>::atom_block(atom_ptr->index, parent->basis, parent->dfbasis);
    int spot = 0;
    nodes[spot]->assign_coef_item(atom_ptr, true);
    while(spot < nodes.size() - 1 and nodes[spot]->estimated_workload > nodes[spot+1]->estimated_workload) {
      std::swap(nodes[spot], nodes[spot+1]);
      ++spot;
    }
  }


}

void
AtomCluster::assign_atom(const boost::shared_ptr<AssignableAtom>& atom) {
  atoms.insert(atom);
  coef_offsets[atom->index] = coefs_size;
  coefs_size += atom->coefs_size;
  dfnsh += atom->dfnshell;
  for(auto&& Xsh : iter_shells_on_center(parent->dfbasis, atom->index, parent->basis)) {
    // Note that atoms must be assigned in order for this to work
    dfbs_shell_map[Xsh.index] = assigned_dfbs_shells.size();
    assigned_dfbs_shells.insert(Xsh.index);
  }
}

void
Assignments::print_detail(std::ostream& o) const
{
  o << indent << "Assignments for exchange computation:" << endl;
  o << incindent;

  o << indent << "Atom clusters:" << endl << incindent;
  for(auto&& cluster : clusters) {
    o << indent << "Cluster " << cluster->index << ": " << endl << incindent;
    o << indent << "Atoms assigned:" << endl << incindent << indent;
    PRINT_AS_GRID(o, atom, cluster->atoms, atom->index, "%5d", 8);
    o << endl << decindent;
    o << indent << "Nodes assigned:" << endl << incindent << indent;
    PRINT_AS_GRID(o, node, cluster->nodes, node->id, "%5d", 8);
    o << endl << decindent;
    o << decindent;
  }
  o << decindent;


  o << indent << "Coefficient memory usage by node:" << endl << incindent << indent;
  PRINT_STRING_AS_GRID(o, node, nodes, (data_size_to_string(node->cluster->coefs_size * sizeof(double))), 8);
  o << endl << decindent;

  o << decindent;

}

