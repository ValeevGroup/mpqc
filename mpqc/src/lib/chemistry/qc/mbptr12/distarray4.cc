//
// distarray4.cc
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/distarray4.h>

using namespace std;
using namespace sc;

namespace {
  // transposes src into dst
  void memcpy_transpose(double* dst, const double* src, int nrow, int ncol) {
    for(int r=0; r<nrow; ++r) {
      double* dst_ptr = dst + r;
      for(int c=0; c<ncol; ++c) {
        *dst_ptr = *src;
        ++src;
        dst_ptr += nrow;
      }
    }
  }
};

DistArray4Dimensions DistArray4Dimensions::default_dim_(-1,-1,-1,-1,-1,DistArray4Storage_XY);
namespace sc {
  bool operator==(const DistArray4Dimensions& A,
                  const DistArray4Dimensions& B) {
    if (&A == &B) return true;
    return A.num_te_types() == B.num_te_types() &&
    A.n1() == B.n1() &&
    A.n2() == B.n2() &&
    A.n3() == B.n3() &&
    A.n4() == B.n4() &&
    A.storage() == B.storage();
  }
}

/*--------------------------------
  DistArray4
 --------------------------------*/
static ClassDesc DistArray4_cd(
  typeid(DistArray4),"DistArray4",2,"virtual public SavableState",
  0, 0, 0);

DistArray4::DistArray4(int num_te_types, int ni, int nj, int nx, int ny,
                       DistArray4Storage storage) :
  num_te_types_(num_te_types), ni_(ni), nj_(nj), nx_(nx), ny_(ny),
  storage_(storage),
  nxy_(nx*ny), blksize_(nxy_*sizeof(double)), blocksize_(blksize_*num_te_types_),
  msg_(MessageGrp::get_default_messagegrp()), active_(false)
{
}

DistArray4::DistArray4(StateIn& si) : SavableState(si),
  msg_(MessageGrp::get_default_messagegrp()), active_(false)
{
  si.get(num_te_types_);
  si.get(ni_);
  si.get(nj_);
  si.get(nx_);
  si.get(ny_);
  int s; si.get(s); storage_ = static_cast<DistArray4Storage>(s);

  nxy_ = nx_ * ny_;
  blksize_ = nxy_*sizeof(double);
  blocksize_ = blksize_*num_te_types_;
}

DistArray4::~DistArray4()
{
}

void DistArray4::save_data_state(StateOut& so)
{
  so.put(num_te_types_);
  so.put(ni_);
  so.put(nj_);
  so.put(nx_);
  so.put(ny_);
  so.put((int)storage_);
}

int
DistArray4::tasks_with_access(vector<int>& twa_map) const
{
  const int nproc = ntasks();

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (has_access(proc)) nproc_with_ints++;

  twa_map.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (has_access(proc)) {
      twa_map[proc] = count;
      count++;
    }
    else
      twa_map[proc] = -1;

  return nproc_with_ints;
}

namespace sc{ namespace detail {

void store_memorygrp(Ref<DistArray4>& acc, Ref<MemoryGrp>& mem, int i_offset,
                     int ni, const size_t blksize_memgrp) {
  // if the accumulator does not accept data from this task, bolt
  if (acc->has_access(mem->me()) == false)
    return;

  // determine over how many tasks the work can be split
  vector<int> writers;
  const int nwriters = acc->tasks_with_access(writers);

  const int me = mem->me();
  const int nproc = mem->n();
  const int num_te_types = acc->num_te_types();
  for (int i=0; i<ni; i++) {
    const int ii = i + i_offset;
    for (int j=0; j<acc->nj(); j++) {
      const int ij = acc->ij_index(i, j);

      // round-robin assignment of blocks among writers
      const int proc_writer = ij % nwriters;
      if (proc_writer != writers[me])
        continue;

      // blocks are distributed in MemoryGrp in round-robin also
      const int proc = ij % nproc;
      const int local_ij_index = ij / nproc;
#if 0
      ExEnv::out0() << indent << "sc::detail::store_memorygrp: me = " << me
                    << " i = " << ii
                    << " j = " << j
                    << " proc = " << proc << endl;
#endif
      if (proc != me) {
        distsize_t moffset = (distsize_t)local_ij_index*blksize_memgrp
            *num_te_types + mem->offset(proc);
        const size_t blksize = acc->blksize();
        for (int te_type = 0; te_type < num_te_types; ++te_type) {
          const double* data = (const double *) mem->obtain_readonly(moffset, blksize);

#if 0
          {
            Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
            RefSCMatrix datamat = localkit->matrix(new SCDimension(acc->nx()),
                                                   new SCDimension(acc->ny()));
            datamat.assign(data);
            std::ostringstream oss;
            oss << "sc::detail::store_memorygrp: i = " << ii << " j = " << j << " te_type = " << te_type;
            datamat.print(oss.str().c_str());
          }
#endif

          acc->store_pair_block(ii, j, te_type, data);
          mem->release_readonly(const_cast<void*>(static_cast<const void*>(data)), moffset, blksize);
          moffset += blksize_memgrp;
        }
      } else {
        const double* data = (const double *) ((size_t)mem->localdata() + blksize_memgrp
            *num_te_types*local_ij_index);
        for (int te_type=0; te_type < num_te_types; te_type++) {

#if 0
          {
            Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
            RefSCMatrix datamat = localkit->matrix(new SCDimension(acc->nx()),
                                                   new SCDimension(acc->ny()));
            datamat.assign(data);
            std::ostringstream oss;
            oss << "sc::detail::store_memorygrp: i = " << ii << " j = " << j << " te_type = " << te_type;
            datamat.print(oss.str().c_str());
          }
#endif

          acc->store_pair_block(ii, j, te_type, data);
          data = (double*) ((size_t) data + blksize_memgrp);
        }
      }
    }
  }
}

void restore_memorygrp(Ref<DistArray4>& acc, Ref<MemoryGrp>& mem, int i_offset,
                       int ni, DistArray4Storage storage, const size_t blksize_memgrp) {
  // if this task cannot get the data from the accumulator, bolt
  if (acc->has_access(mem->me()) == false)
    return;

  const bool transpose_xy = (storage != acc->storage());
  int nrow=-1, ncol=-1;
  if (transpose_xy) {
    // number of row indices in blocks as stored in acc
    nrow = (acc->storage() == DistArray4Storage_XY) ? acc->nx() : acc->ny();
    // number of column indices in blocks as stored in acc
    ncol = (acc->storage() == DistArray4Storage_XY) ? acc->ny() : acc->nx();
  }

  // determine over how many tasks the work can be split
  vector<int> readers;
  const int nreaders = acc->tasks_with_access(readers);

  const int me = mem->me();
  const int nproc = mem->n();
  const int num_te_types = acc->num_te_types();
  for (int i=0; i<ni; i++) {
    const int ii = i + i_offset;
    for (int j=0; j<acc->nj(); j++) {
      const int ij = acc->ij_index(i, j);

      // round-robin assignment of blocks among readers
      const int proc_reader = ij % nreaders;
      if (proc_reader != readers[me])
        continue;

      // blocks are distributed in MemoryGrp in round-robin also
      const int proc = ij % nproc;
      const int local_ij_index = ij / nproc;
      if (proc != me) {
        distsize_t moffset = (distsize_t)local_ij_index*blksize_memgrp
            *num_te_types + mem->offset(proc);
        const size_t blksize = acc->blksize();
        for (int te_type = 0; te_type < num_te_types; ++te_type) {
          const double* data = acc->retrieve_pair_block(ii, j, te_type);
          double* buffer = (double *) mem->obtain_writeonly(moffset, blksize);

          if (!transpose_xy)
            ::memcpy((void*)buffer, (const void*)data, blksize);
          else {
            memcpy_transpose(buffer, data, nrow, ncol);
          }

          mem->release_writeonly(const_cast<void*>(static_cast<const void*>(buffer)), moffset, blksize);
          moffset += blksize_memgrp;
        }
      } else {
        double* buffer = (double *) ((size_t)mem->localdata() + blksize_memgrp
            *num_te_types*local_ij_index);
        const size_t blksize = acc->blksize();
        for (int te_type=0; te_type < num_te_types; te_type++) {
          const double* data = acc->retrieve_pair_block(ii, j, te_type);

          if (!transpose_xy)
            ::memcpy((void*)buffer, (const void*)data, blksize);
          else {
            memcpy_transpose(buffer, data, nrow, ncol);
          }

          acc->release_pair_block(ii, j, te_type);
          buffer = (double*) ((size_t) buffer + blksize_memgrp);
        }
      }
    }
  }
}

}} // end of namespace sc::detail

namespace sc {

  /** Algorithm:
   *  given src = (ij|xy), clone to create (ix|jy), t is the number of te types
   *  tile x such that the entire jy block can be held (X*nj*ny <= memory, where Y is the tilesize)
   *  loop over all i, t
   *    loop over x tiles
   *      assign i,t,xtile to a worker
   *      loop over all j
   *        read xtile,y block
   *        sort to x,j,y
   *      end j loop
   *      loop over x in this tile
   *        write jy
   *      end x loop
   *    end x tile loop
   *  end i, t loop
   */
  Ref<DistArray4> permute23(const Ref<DistArray4>& src,
                            size_t available_memory) {

    const int nt = src->num_te_types();
    const int ni = src->ni();
    const int nj = src->nj();
    const int nx = src->nx();
    const int ny = src->ny();
    const int njy = nj * ny;
    DistArray4Dimensions result_dims(nt, ni, nx, nj, ny, DistArray4Storage_XY);
    Ref<DistArray4> result = src->clone(result_dims);
    result->activate();

    // determine the size and number of x tiles
    const size_t jy_blksize = njy * sizeof(double);
    int tilesize = (available_memory + jy_blksize - 1) / jy_blksize;
    if (tilesize > nx)  tilesize = nx;
    const int ntiles = (nx + tilesize - 1) / tilesize;
    tilesize = (nx + ntiles - 1) / ntiles;  // recompute the tile size for a given number of tiles
                                            // to make the tiles more even

    // allocate memory for the work buffers
    double* xjy_buf = new double[tilesize * njy];
    double* xy_buf = new double[tilesize * ny];

    // determine how many workers we have
    std::vector<int> worker_id;
    const int nworkers = src->tasks_with_access(worker_id);
    const int me = MessageGrp::get_default_messagegrp()->me();

    // sort!
    int task_id = 0;
    for(int i=0; i<ni; ++i) {
      for(int t=0; t<nt; ++t) {

        int xoffset = 0;
        for(int tile=0; tile<ntiles; ++tile, xoffset+=tilesize, ++task_id){

          // round-robin task allocation
          if(task_id%nworkers != worker_id[me])
            continue;

          int xfence = xoffset+tilesize;
          if (xfence > nx) xfence = nx;
          const int this_tilesize = xfence - xoffset;
          for(int j=0; j<nj; ++j) {
            src->retrieve_pair_subblock(i, j, t,
                                        xoffset, xfence, 0, ny,
                                        xy_buf);

            double* xj_ptr = xjy_buf + j*ny;
            const double* x_ptr = xy_buf;
            // copy each y row to its location in xjy_buf
            for(int x=0; x<this_tilesize; ++x, x_ptr+=ny, xj_ptr+=njy) {
              std::copy(x_ptr, x_ptr+ny, xj_ptr);
            }
          }

          // write i,x blocks to the destination
          const double* x_ptr = xjy_buf;
          for(int x=0; x<this_tilesize; ++x, x_ptr+=njy) {
            result->store_pair_block(i, x+xoffset, t, x_ptr);
          }

        }

      }
    }

    delete[] xjy_buf;
    delete[] xy_buf;

    return result;
  }

} // end of namespace sc

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
