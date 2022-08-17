
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_parallel_h
#define _chemistry_qc_lmp2_parallel_h

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#include <memory>
#include <algorithm>
#include <random>
#include <mpqc_config.h>
#include <util/misc/scint.h>
#include <util/misc/regtime.h>

#include <chemistry/qc/lmp2/collective.h>

#define NEW_DIST_DATA

namespace sc {

namespace sma2 {

template <int N>
void
Array<N>::parallel_accumulate(const sc::Ref<sc::MessageGrp> &grp) {
  sc::Timer tim;
  tim.enter("accum");
  //tim_enter("sum");
  //grp->sum(data_->data(), data_->ndata());
  //tim_exit("sum");
#ifndef HAVE_MPI
  tim.enter("sum0");
  grp->sum(data_->data(), data_->ndata(), 0, 0);
  tim.exit("sum0");
  tim.enter("bcast");
  grp->bcast(data_->data(), data_->ndata());
  tim.exit("bcast");
#else
  allsum(data_->data(), data_->ndata());
#endif
  tim.enter("bounds");
  compute_bounds();
  tim.exit("bounds");
  tim.exit("accum");
}

template <int N>
void
Array<N>::parallel_union(const sc::Ref<sc::MessageGrp> &msg)
{
  sc::Timer tim;
  int me = msg->me();
  int nproc = msg->n();

  // fetch the block map entries from all the other nodes
  tim.enter("parallel_union");
  int nblock_local = n_block();
  int nblock_max = n_block();
  msg->max(nblock_max);
  int *local_blocks = new int[1 + nblock_local*N];
  int *remote_blocks = new int[1 + nblock_max*N];
  int iarray = 0;
  local_blocks[iarray++] = nblock_local;
  for (typename blockmap_t::iterator i = blocks_.begin();
       i != blocks_.end();
       i++) {
      for (int j=0; j<N; j++) {
          local_blocks[iarray++] = i->first.block(j);
        }
    }
  int source = (me + nproc - 1)%nproc; //  me - 1 mod nproc
  for (int next_node = (me+1)%nproc;
       next_node != me;
       next_node = (next_node+1)%nproc) {
#ifdef HAVE_MPI
      MPI_Request send_req, recv_req;
      MPI_Irecv(remote_blocks,1+nblock_max*N,MPI_INT,source,
                100+source,MPI_COMM_WORLD,&recv_req);
      MPI_Isend(local_blocks,1+nblock_local*N,MPI_INT,next_node,
                100+me,MPI_COMM_WORLD,&send_req);
      MPI_Status status;
      MPI_Wait(&recv_req,&status); // wait for the recv to complete
      int iarray = 0;
      int nblock = remote_blocks[iarray++];
      for (int i=0; i<nblock; i++) {
          BlockInfo<N> bi;
          for (int j=0; j<N; j++) bi.block(j) = remote_blocks[iarray++];
          add_unallocated_block(bi);
        }
      MPI_Wait(&send_req,&status); // wait for the send to complete
      source = (source + nproc - 1)%nproc; // source - 1 mod nproc
#else
      throw std::runtime_error("Array<N,C>::replicated_from_distributed requires MPI");
#endif
    }
  delete[] local_blocks;
  delete[] remote_blocks;
  tim.exit("parallel_union");
}

/** \brief Provides information about how blocks are distributed onto processes.
 */
template <int N>
class BlockDistrib {
  public:
    /// Given a block, returns the node on which it resides.
    virtual int block_to_node(const BlockInfo<N> &b) const = 0;
    virtual ~BlockDistrib() {}
};

/** \brief Distributes pairs of indices among the processes.
 */
class PairMapping: public sc::RefCount {
    std::map<std::pair<int,int>,int> nodemap_;
    int nproc_;
    int me_;
  public:
    // Distributes two indices among the nodes.  The passed
    // indexset must be identical on all nodes.
    PairMapping(const sc::Ref<sc::MessageGrp> &grp,
                const std::set<std::pair<int,int> > &indexset,
                bool randomize = false) {
      nproc_ = grp->n();
      me_= grp->me();
      std::vector<int> nodemap(indexset.size());
      for (int i=0; i<nodemap.size(); i++) nodemap[i] = i%nproc_;
      if (randomize) {
          std::shuffle(nodemap.begin(), nodemap.end(), std::default_random_engine(0));
        }
      int node=0;
      for (std::set<std::pair<int,int> >::const_iterator i = indexset.begin();
           i != indexset.end();
           i++) {
          nodemap_[*i] = nodemap[node++];
        }
    }
    // Distributes indices i0 and i1 among the nodes.  The passed
    // indexmap is identical on all nodes.
    PairMapping(const sc::Ref<sc::MessageGrp> &grp,
                const std::map<int,std::set<int> > &indexmap) {
      nproc_ = grp->n();
      me_= grp->me();
      int node=0;
      for (std::map<int,std::set<int> >::const_iterator i = indexmap.begin();
           i != indexmap.end();
           i++) {
          for (std::set<int>::const_iterator j = i->second.begin();
               j != i->second.end();
               j++) {
              nodemap_[std::make_pair(i->first,*j)] = node++%nproc_;
            }
        }
    }
    int pair_to_node(int i, int j) const {
      std::map<std::pair<int,int>,int>::const_iterator
          found = nodemap_.find(std::make_pair(i,j));
      if (found == nodemap_.end()) {
          throw std::runtime_error("PairMapping: requested non-existent block");
        }
      return found->second;
    }
    int me() const { return me_; }
    int nproc() const { return nproc_; }
    void local_pairs(std::set<std::pair<int,int> > &pairs) const {
      pairs.clear();
      for (std::map<std::pair<int,int>,int>::const_iterator i=nodemap_.begin();
           i != nodemap_.end();
           i++) {
          if (i->second == me_) pairs.insert(i->first);
        }
    }
    void print() const {
      for (std::map<std::pair<int,int>,int>::const_iterator i=nodemap_.begin();
           i != nodemap_.end();
           i++) {
          sc::ExEnv::out0() << i->first.first
                            << " " << i->first.second
                            << " -> " << i->second
                            << std::endl;
        }
    }
};

/** \brief An implementation of BlockDistrib using PairMapping.
 */
template <int N>
class PairBlockDistrib: public BlockDistrib<N> {
    int i0_, i1_;
    sc::Ref<PairMapping> mapping_;
  public:
    PairBlockDistrib(const sc::Ref<sc::MessageGrp> &grp,
                     int i0, int i1,
                     const std::set<std::pair<int,int> > &indexset,
                     bool randomize = false) {
      i0_ = i0; i1_ = i1;
      mapping_ = new PairMapping(grp,indexset,randomize);
    }
    PairBlockDistrib(const sc::Ref<sc::MessageGrp> &grp,
                     int i0, int i1,
                     const std::map<int,std::set<int> > &indexmap) {
      i0_ = i0; i1_ = i1;
      mapping_ = new PairMapping(grp,indexmap);
    }
    // Distributes indices i0 and i1 among the nodes.
    // The PairMapping object describes the distribution.
    PairBlockDistrib(int i0, int i1,
                     const sc::Ref<PairMapping> &mapping) {
      i0_ = i0; i1_ = i1;
      mapping_ = mapping;
    }
    int block_to_node(const BlockInfo<N> &b) const {
      return mapping_->pair_to_node(b.block(i0_),b.block(i1_));
    }
    void local_pairs(std::set<std::pair<int,int> > &pairs) const {
      mapping_->local_pairs(pairs);
    }
};

/** \brief Distribute blocks round-robin among processes using one or
    more index values.
 */
template <int N>
class CompleteBlockDistrib: public BlockDistrib<N> {
    int nindex_;
    int indices_[N];
    sc::sc_uint64_t size_[N];
    int nindex_i_[N];
    int nproc_;
    int me_;
    void init_sizes(const Array<N> &a) {
      size_[0] = 1;
      for (int i=1; i<N; i++) {
          size_[i] = size_[i-1] * a.index(i).nindex();
        }
      for (int i=0; i<N; i++) {
          nindex_i_[i] = a.index(i).nindex();
        }
    }
  public:
    CompleteBlockDistrib(const Array<N> &a, const sc::Ref<sc::MessageGrp> &grp,
                         int i0): nindex_(1) {
      me_ = grp->me();
      indices_[0]=i0;
      nproc_ = grp->n();
      init_sizes(a);
    }
    CompleteBlockDistrib(const Array<N> &a, const sc::Ref<sc::MessageGrp> &grp,
                         int i0, int i1): nindex_(2) {
      me_ = grp->me();
      indices_[0]=i0; indices_[1]=i1;
      nproc_ = grp->n();
      init_sizes(a);
    }
    CompleteBlockDistrib(const Array<N> &a, const sc::Ref<sc::MessageGrp> &grp,
                         int i0, int i1, int i2): nindex_(3) {
      me_ = grp->me();
      indices_[0]=i0; indices_[1]=i1; indices_[2]=i2;
      nproc_ = grp->n();
      init_sizes(a);
    }
    int block_to_node(const BlockInfo<N> &b) const {
      sc::sc_uint64_t loc = 0;
      for (int i=0; i<nindex_; i++) {
          loc += size_[i] * b.block(indices_[i]);
        }
      return loc%nproc_;
    }
    // Only valid for a single index!
    int nlocalindex(int node = -1) {
      if (node == -1) node = me_;
      int r = nindex_i_[0] / nproc_;
      if (me_ < nindex_i_[0] % nproc_) r++;
      return r;
    }
    // Only valid for a single index!
    int localindex(int i, int node = -1) {
      if (node == -1) node = me_;
      return me_ + nproc_*i;
    }
};

template <int N>
void
Array<N>::replicated_from_distributed(const sc::Ref<sc::MessageGrp>& msg,
                                      const Array<N> &d)
{
  sc::Timer tim;
  tim.enter("repl_from_dist");

  int me = msg->me();
  int nproc = msg->n();
//   if (nproc == 1) {
//       operator = (d);
//       return;
//     }

  // initialize the indices, etc, in the array
  // and copy the blockmap entries from the local node
  init_blocks(d, 0.0); // using tol=0.0 causes all blocks to be used
#ifdef USE_BOUND
  bound_ = d.bound_;
  msg->max(bound_);
#endif
  tolerance_ = d.tolerance_;

  parallel_union(msg);

  // allocate and initialize storage
  zero();

  // copy the data from the local node
  tim.enter("local");
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int ndat = d.block_size(i->first);
      typename blockmap_t::iterator loc = blocks_.find(i->first);
      // must make sure block was found, since init_blocks may
      // throw out blocks with a bound below the tolerance
      if (loc != blocks_.end()) {
          memcpy(loc->second, i->second, ndat*sizeof(double));
        }
    }
  tim.exit("local");

  // accumulate the data on all nodes
  parallel_accumulate(msg);

  tim.exit("repl_from_dist");
}

template <int N>
void
Array<N>::distributed_from_distributed(const sc::Ref<sc::MessageGrp>& msg,
                                       const BlockDistrib<N> &distrib,
                                       Array<N> &d,
                                       bool clear_source_array,
                                       bool ignore_block_distrib_throws)
{
  int me = msg->me();
  int nproc = msg->n();

  // If MPI exists, then the MPI routines will be used on even a single
  // processor.  This is done for testing purposes.
#ifndef HAVE_MPI
  if (nproc == 1) {
      assign_all(d);
      if (clear_source_array) d.clear();
      return;
    }
  else {
      throw std::runtime_error("Array<N,C>::distributed_from_distributed requires MPI");
    }
#else // HAVE_MPI
  // clj debug
  //d.print(msg,true,std::cout);
  // clj end debug

  // initialize the array
  clear();
  set_indices(d.indices());

#if 1
  // compute number of blocks and amount of data moving from this node to
  // each other node
  int *nblock_send = new int[nproc];
  int *nblock_recv = new int[nproc];
  for (int i=0; i<nproc; i++) nblock_send[i] = 0;
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int node;
      if (ignore_block_distrib_throws) {
          try {
              node = distrib.block_to_node(i->first);
            }
          catch (...) {
              continue;
            }
        }
      else {
          node = distrib.block_to_node(i->first);
        }
      nblock_send[node]++;
    }

  // Exchange the size of buffers that will be sent from each node.
  MPI_Alltoall(nblock_send, 1, MPI_INT,
               nblock_recv, 1, MPI_INT,
               MPI_COMM_WORLD);

  // clj debug
//   for (int i=0; i<nproc; i++) {
//       std::cout << me << "nblock_send[" << i << "] = " << nblock_send[i]
//                 << std::endl;
//     }
//   for (int i=0; i<nproc; i++) {
//       std::cout << me << "nblock_recv[" << i << "] = " << nblock_recv[i]
//                 << std::endl;
//     }
  // clj end debug

  // Find out how much data we'll send and recv and what the offsets are.
  int *block_send_offset = new int[nproc];
  int *block_recv_offset = new int[nproc];
  int *block_send_size = new int[nproc];
  int *block_recv_size = new int[nproc];
  for (int i=0; i<nproc; i++) {
      block_send_size[i] = N*nblock_send[i];
      block_recv_size[i] = N*nblock_recv[i];
      if (i==0) {
          block_send_offset[i] = 0;
          block_recv_offset[i] = 0;
        }
      else {
          block_send_offset[i] = block_send_size[i-1]
                               + block_send_offset[i-1];
          block_recv_offset[i] = block_recv_size[i-1]
                               + block_recv_offset[i-1];
        }
    }
  int block_send_size_total = block_send_size[nproc-1]
                            + block_send_offset[nproc-1];
  int block_recv_size_total = block_recv_size[nproc-1]
                            + block_recv_offset[nproc-1];

  // Allocate the blockinfo buffers and copy data into them.  The ordering
  // of data is nblock, blockinfo[0], ...,
  // blockinfo[n-1].
  int *block_send_buffers = new int[block_send_size_total];
  std::vector<int> current_send_offset(nproc);
  std::fill(current_send_offset.begin(), current_send_offset.end(), 0);
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int node;
      if (ignore_block_distrib_throws) {
          try {
              node = distrib.block_to_node(i->first);
            }
          catch (...) {
              continue;
            }
        }
      else {
          node = distrib.block_to_node(i->first);
        }
      for (int j=0; j<N; j++)
          block_send_buffers[block_send_offset[node]
                             + current_send_offset[node]++]
              = i->first.block(j);
    }

  // Exchange the blockinfo buffers.
  int *block_recv_buffers = new int[block_recv_size_total];
  custom_alltoallv(block_send_buffers, block_send_size,
                block_send_offset, MPI_INT,
                block_recv_buffers, block_recv_size,
                block_recv_offset, MPI_INT,
                MPI_COMM_WORLD);

  // Add the blockinfos to this Array object and allocate the data.
  int irecv = 0;
  long ndata_recv_total = 0;
  int *data_recv_offset = new int[nproc];
  int *data_recv_size = new int[nproc];
  for (int i=0; i<nproc; i++) {
      int nblock = nblock_recv[i];
      long ndata = 0;
      for (int j=0; j<nblock; j++) {
          BlockInfo<N> bi;
          for (int k=0; k<N; k++) {
              bi.block(k) = block_recv_buffers[irecv++];
            }
          ndata += block_size(bi);
          add_unallocated_block(bi);
        }
      ndata_recv_total += ndata;
      data_recv_size[i] = ndata;
      if (i == 0)
          data_recv_offset[i] = 0;
      else
          data_recv_offset[i] = data_recv_offset[i-1] + data_recv_size[i-1];
    }

  // Allocate the buffers for the exchange of data and place the data into
  // the send buffer.
  double *data_send_buffers = new double[d.n_element_allocated()];
  int *data_send_offset = new int[nproc];
  int *data_send_size = new int[nproc];
  int isend = 0;
  double *current_send_buffer = data_send_buffers;
  for (int i=0; i<nproc; i++) {
      int nblock = nblock_send[i];
      data_send_size[i] = 0;
      for (int j=0; j<nblock; j++) {
          BlockInfo<N> bi;
          for (int k=0; k<N; k++) {
              bi.block(k) = block_send_buffers[isend++];
            }
          int ndata_in_block = block_size(bi);
          memcpy(current_send_buffer, d.blocks_.find(bi)->second,
                 ndata_in_block*sizeof(double));
          current_send_buffer += ndata_in_block;
          data_send_size[i] += ndata_in_block;
        }
      if (i == 0)
          data_send_offset[i] = 0;
      else
          data_send_offset[i] = data_send_offset[i-1] + data_send_size[i-1];
    }
  delete[] block_send_buffers;
  delete[] block_send_offset;
  delete[] block_send_size;
  if (clear_source_array) d.clear();

  // Exchange the data
  double *data_recv_buffers = new double[ndata_recv_total];
  custom_alltoallv(data_send_buffers, data_send_size,
                data_send_offset, MPI_DOUBLE,
                data_recv_buffers, data_recv_size,
                data_recv_offset, MPI_DOUBLE,
                MPI_COMM_WORLD);
  delete[] data_send_buffers;
  delete[] data_send_offset;
  delete[] data_send_size;

  // Update the array's data.
  allocate_blocks();
  irecv = 0;
  double *current_recv_buffer = data_recv_buffers;
  for (int i=0; i<nproc; i++) {
      int nblock = nblock_recv[i];
      for (int j=0; j<nblock; j++) {
          BlockInfo<N> bi;
          for (int k=0; k<N; k++) {
              bi.block(k) = block_recv_buffers[irecv++];
            }
          typename blockmap_t::iterator iblock = blocks_.find(bi);
          double *data = iblock->second;
          int ndata_in_block = block_size(bi);
          memcpy(data, current_recv_buffer, ndata_in_block*sizeof(double));
          current_recv_buffer += ndata_in_block;
        }
    }

  delete[] block_recv_buffers;
  delete[] block_recv_offset;
  delete[] block_recv_size;

  delete[] data_recv_buffers;
  delete[] data_recv_offset;
  delete[] data_recv_size;

  delete[] nblock_send;
  delete[] nblock_recv;

#else

  /////////////////////////////////////////////////////////////////////
  // fetch my block map entries from all the other nodes
  // and send d's contributions to all the other nodes

  // compute number of blocks moving from this node to each other node
  std::vector<int> nblock(nproc);
  std::fill(nblock.begin(), nblock.end(), 0);
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int node = distrib.block_to_node(i->first);
      nblock[node]++;
    }

  // allocate buffers for each node
  std::vector<int*> buffer(nproc);
  for (int i=0; i<nproc; i++) {
      if (i != me) {
          buffer[i] = new int[1 + N*nblock[i]];
          buffer[i][0] = nblock[i];
        }
      else {
          buffer[i] = 0;
        }
    }

  // pack the buffer for each node with blockinfo data
  // local blocks are directly inserted
  std::vector<int> ndata(nproc);
  std::fill(ndata.begin(), ndata.end(), 1); // fill w/1 since all contain nblocks 1st
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int node = distrib.block_to_node(i->first);
      int &counter = ndata[node];
      for (int j=0; j<N; j++) {
          if (node == me) {
              add_unallocated_block(i->first);
            }
          else {
              buffer[node][counter++] = i->first.block(j);
            }
        }
    }

  // send the blockinfo data from this node to each other node
  int max_nblock = *std::max_element(nblock.begin(), nblock.end());
  msg->max(max_nblock);
  int *remote_blocks = new int[1 + max_nblock*N];
  int source = (me + nproc - 1)%nproc; //  me - 1 mod nproc
  for (int next_node = (me+1)%nproc;
       next_node != me;
       next_node = (next_node+1)%nproc) {
#ifdef HAVE_MPI
      MPI_Request send_req, recv_req;
      MPI_Irecv(remote_blocks,1+max_nblock*N,MPI_INT,source,
                100+source,MPI_COMM_WORLD,&recv_req);
      MPI_Isend(buffer[next_node],1+nblock[next_node]*N,MPI_INT,next_node,
                100+me,MPI_COMM_WORLD,&send_req);
      MPI_Status status;
      MPI_Wait(&recv_req,&status); // wait for the recv to complete
      int iarray = 0;
      int nblock = remote_blocks[iarray++];
      for (int i=0; i<nblock; i++) {
          BlockInfo<N> bi;
          for (int j=0; j<N; j++) bi.block(j) = remote_blocks[iarray++];
          add_unallocated_block(bi);
        }
      MPI_Wait(&send_req,&status); // wait for the send to complete
      source = (source + nproc - 1)%nproc; // source - 1 mod nproc
#else
      throw std::runtime_error("Array<N,C>::distributed_from_distributed requires MPI");
#endif
    }
  delete[] remote_blocks;
  for (int i=0; i<nproc; i++) delete[] buffer[i];

  /////////////////////////////////////////////////////////////////////
  // allocate storage
  allocate_blocks();
  // find the largest packet size and allocate storage
  int maxsz = 0;
#ifdef HAVE_MPI
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end();
       i++) {
      int sz = block_size(i->first);
      if (sz > maxsz) maxsz = sz;
    }
  msg->max(maxsz);
  int maxalloc;
  MPI_Pack_size(maxsz, MPI_DOUBLE, MPI_COMM_WORLD, &maxalloc);
  for (int i=0; i<N; i++) {
      int isz;
      MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &isz);
      maxalloc += isz;
    }
  char *outgoing_buffer = new char[maxalloc];
#endif

  /////////////////////////////////////////////////////////////////////
  // send/recv data
  int n_to_be_received = blocks_.size() - nblock[me];
#ifdef HAVE_MPI
  const int nrecv_req = 20;
  std::vector<char *> incoming_buffer(nrecv_req);
  for (int i=0; i<nrecv_req; i++) {
      incoming_buffer[i] = new char[maxalloc];
    }
  MPI_Request recv_req[nrecv_req];
  for (int i=0; i<nrecv_req; i++) {
      if (i<n_to_be_received) {
          MPI_Irecv(incoming_buffer[i], maxalloc, MPI_BYTE, MPI_ANY_SOURCE,
                    99, MPI_COMM_WORLD, &recv_req[i]);
        }
      else {
          recv_req[i] = MPI_REQUEST_NULL;
        }
    }
#endif
  for (typename blockmap_t::const_iterator i = d.blocks_.begin();
       i != d.blocks_.end() || n_to_be_received > 0;) {
#ifdef HAVE_MPI
      MPI_Request send_req;
#endif
      int node;
//       std::cout << sc::scprintf("%d end = %d nleft = %d\n",
//                                 me, int(i == d.blocks_.end()), n_to_be_received)
//                 << std::flush;
      bool need_wait_for_send = false;
      if (i != d.blocks_.end()) {
          node = distrib.block_to_node(i->first);
          int ndat = d.block_size(i->first);
          if (node == me) {
//               std::cout << sc::scprintf("  %d processing local block\n",me)
//                         << std::flush;
              typename blockmap_t::iterator loc = blocks_.find(i->first);
              memcpy(loc->second, i->second, ndat*sizeof(double));
            }
          else {
//               std::cout << sc::scprintf("  %d sending block to %d\n",me,node)
//                         << std::flush;
              // pack up block and send it off
#ifdef HAVE_MPI
              int outgoing_loc = 0;
              for (int j=0; j<N; j++) {
                  int jval = i->first.block(j);
                  MPI_Pack(&jval, 1, MPI_INT, outgoing_buffer, maxalloc,
                           &outgoing_loc, MPI_COMM_WORLD);
                }
              MPI_Pack(i->second, ndat, MPI_DOUBLE, outgoing_buffer, maxalloc,
                       &outgoing_loc, MPI_COMM_WORLD);
              MPI_Isend(outgoing_buffer,outgoing_loc,MPI_BYTE,node,
                        99,MPI_COMM_WORLD,&send_req);
              need_wait_for_send = true;
#else
              throw std::runtime_error("LMP2:: cannot send remote block w/o MPI");
#endif
            }
          i++;
        }
#ifdef HAVE_MPI
      // See if any data is ready to be recved
      int nrecvd = 0;
      for (int irecv=0; irecv<nrecv_req; irecv++) {
          if (recv_req[irecv] == MPI_REQUEST_NULL) continue;
          int flag = 1;
          MPI_Status status;
          MPI_Test(&recv_req[irecv],&flag,&status);
          if (flag) {
//               std::cout << sc::scprintf("  %d got incoming\n",me)
//                         << std::flush;
              BlockInfo<N> bi;
              int position = 0;
              for (int j=0; j<N; j++) {
                  int index;
                  MPI_Unpack(incoming_buffer[irecv], maxalloc, &position, &index,
                             1, MPI_INT, MPI_COMM_WORLD);
                  bi.block(j) = index;
                }
              typename blockmap_t::iterator loc = blocks_.find(bi);
              int block_sz = block_size(bi);
              MPI_Unpack(incoming_buffer[irecv], maxalloc, &position, loc->second,
                         block_sz, MPI_DOUBLE, MPI_COMM_WORLD);
              n_to_be_received--;
              // only repost the receive if there are not enough recv requests remaining
              if (n_to_be_received >= nrecv_req) {
                  MPI_Irecv(incoming_buffer[irecv], maxalloc, MPI_BYTE, MPI_ANY_SOURCE,
                            99, MPI_COMM_WORLD, &recv_req[irecv]);
                }
            }
        }

      // Wait until data is sent
      if (need_wait_for_send) {
          MPI_Status status;
          MPI_Wait(&send_req,&status);
        }
#endif
    }
#ifdef HAVE_MPI
  for (int i=0; i<nrecv_req; i++) {
      delete[] incoming_buffer[i];
    }
  delete[] outgoing_buffer;
#endif
  if (clear_source_array) d.clear();
#endif
#endif // HAVE_MPI

  compute_bounds();

  // clj debug
  //print(msg,true,std::cout);
  // clj end debug
}

}

}

#endif
