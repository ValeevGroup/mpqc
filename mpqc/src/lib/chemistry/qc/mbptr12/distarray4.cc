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
#include <util/misc/consumableresources.h>
#include <util/misc/regtime.h>
#include <util/class/scexception.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/print.h>

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

  void contract34(const Ref<DistArray4>& braket,
                  double scale,
                  const Ref<DistArray4>& bra,
                  unsigned int intsetidx_bra,
                  const Ref<DistArray4>& ket,
                  unsigned int intsetidx_ket,
                  int debug) {

    bra->activate();
    ket->activate();
    braket->activate();

    const unsigned int nb1 = bra->ni();
    const unsigned int nb2 = bra->nj();
    const unsigned int nk1 = ket->ni();
    const unsigned int nk2 = ket->nj();
    const unsigned int n1 =  bra->nx();
    const unsigned int n2 =  bra->ny();
    assert(n1 == ket->nx());
    assert(n2 == ket->ny());
    assert(braket->ni() == nb1);
    assert(braket->nj() == nb2);
    assert(braket->nx() == nk1);
    assert(braket->ny() == nk2);

    // Using spinorbital iterators means I don't take into account perm symmetry
    // More efficient algorithm will require generic code
    SpinMOPairIter iterbra(nb1, nb2, false);
    SpinMOPairIter iterket(nk1, nk2, false);
    SpinMOPairIter iterint( n1,  n2, false);
    // size of one block of <space1_bra space2_bra|
    const unsigned int nbra = iterbra.nij();
    // size of one block of <space1_ket space2_ket|
    const unsigned int nket = iterket.nij();

    //
    // data will be accessed in tiles
    // determine tiling here
    //
    // total number of (bra1 bra2| index combinations
    const size_t nij_bra = nbra;
    // total number of (ket1 ket2| index combinations
    const size_t nij_ket = nket;
    // size of each integral block |int1 int2)
    const unsigned int blksize_int_sq = n1 * n2;
    const unsigned int blksize_int = blksize_int_sq;
    //
    // maximum tile size is determined by the available memory
    const size_t memory_available = ConsumableResources::get_default_instance()->memory();
    const size_t max_tile_size = memory_available / (2 * blksize_int * sizeof(double));
    if (max_tile_size == 0) {
      throw AlgorithmException("not enough memory for a single tile, increase memory", __FILE__, __LINE__);
    }
    // try tiling nij_bra and nij_ket into nproc tiles
    const int nproc = bra->msg()->n();
    size_t try_tile_size_bra = (nij_bra + nproc - 1) / nproc;
    size_t try_tile_size_ket = (nij_ket + nproc - 1) / nproc;
    try_tile_size_bra = std::min(try_tile_size_bra, max_tile_size);
    try_tile_size_ket = std::min(try_tile_size_ket, max_tile_size);
    const size_t ntiles_bra = (nij_bra + try_tile_size_bra - 1) / try_tile_size_bra;
    const size_t ntiles_ket = (nij_ket + try_tile_size_ket - 1) / try_tile_size_ket;
    const size_t tile_size_bra = (nij_bra + ntiles_bra - 1) / ntiles_bra;
    const size_t tile_size_ket = (nij_ket + ntiles_ket - 1) / ntiles_ket;

    // scratch buffers to hold tiles and the result of the contraction
    double* T_bra = new double[tile_size_bra * blksize_int];
    double* T_ket = new double[tile_size_ket * blksize_int];
    double* T_result = new double[tile_size_bra * tile_size_ket];

    // split work over tasks which have access to integrals
    // WARNING: assuming same accessibility for both bra and ket transforms
    std::vector<int> proc_with_ints;
    const int nproc_with_ints = bra->tasks_with_access(proc_with_ints);
    const int me = bra->msg()->me();

    if (bra->has_access(me)) {

      size_t task_count = 0;

      // these will keep track of each tile's sets of i and j values
      typedef detail::triple<unsigned int, unsigned int, unsigned int> uint3;
      std::vector<uint3> bra_ij;
      iterbra.start();
      // loop over bra tiles for this set
      size_t tbra_offset = 0;
      for (size_t tbra = 0; tbra < ntiles_bra; ++tbra, tbra_offset += tile_size_bra) {
        double* bra_tile = 0;   // zero indicates it needs to be loaded

        // these will keep track of each tile's sets of i and j values
        std::vector<uint3> ket_ij;
        iterket.start();
        // loop over ket tiles for this set
        size_t tket_offset = 0;
        for (size_t tket = 0; tket < ntiles_ket; ++tket, tket_offset += tile_size_ket, ++task_count) {
          double* ket_tile = 0;   // zero indicates it needs to be loaded

          // distribute tasks by round-robin
          const int task_proc = task_count % nproc_with_ints;
          if (task_proc != proc_with_ints[me])
            continue;

          // has the bra tile been loaded?
          if (bra_tile == 0) {
            bra_tile = T_bra;
            for (size_t i=0; iterbra && i<tile_size_bra; ++i, iterbra.next()) {
              const uint3 ijt(iterbra.i(), iterbra.j(), iterbra.ij());
              bra_ij.push_back(ijt);

              double* blk_ptr = bra_tile + i*blksize_int;
              // where to read the integrals? if not antisymmetrizing read directly to T_bra, else read to scratch
              // buffer, then antisymmetrize to T_bra
              double* blk_ptr_read = blk_ptr;
              Timer tim_intsretrieve("MO ints retrieve");
              const double* blk_cptr = bra->retrieve_pair_block(ijt.i0_,
                                                                   ijt.i1_,
                                                                   intsetidx_bra,
                                                                   blk_ptr_read);
              tim_intsretrieve.exit();

              if (debug >= DefaultPrintThresholds::allO2N2) {
                ExEnv::outn() << indent << "task " << me
                    << ": obtained ij blocks" << std::endl;
                ExEnv::outn() << indent
                              << "i = " << ijt.i0_
                              << " j = " << ijt.i1_ << std::endl;

                RefSCMatrix blk_scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(n1),
                                                                                 new SCDimension(n2));
                blk_scmat.assign(blk_cptr);
                blk_scmat.print("ij block");
              }

            }
          }

          // has the ket tile been loaded?
          if (ket_tile == 0) {
            ket_tile = T_ket;
            for (size_t i=0; iterket && i<tile_size_ket; ++i, iterket.next()) {
              const uint3 ijt(iterket.i(), iterket.j(), iterket.ij());
              ket_ij.push_back(ijt);

              double* blk_ptr = ket_tile + i*blksize_int;
              // where to read the integrals? if not antisymmetrizing read directly to T_ket, else read to scratch
              // buffer, then antisymmetrize to T_ket
              double* blk_ptr_read = blk_ptr;
              Timer tim_intsretrieve("MO ints retrieve");
              const double* blk_cptr = ket->retrieve_pair_block(ijt.i0_,
                                                                   ijt.i1_,
                                                                   intsetidx_ket,
                                                                   blk_ptr_read);
              tim_intsretrieve.exit();

              if (debug >= DefaultPrintThresholds::allO2N2) {
                ExEnv::outn() << indent << "task " << me
                    << ": obtained kl blocks" << std::endl;
                ExEnv::outn() << indent
                              << "k = " << ijt.i0_
                              << " l = " << ijt.i1_ << std::endl;

                RefSCMatrix blk_scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(n1),
                                                                                 new SCDimension(n2));
                blk_scmat.assign(blk_cptr);
                blk_scmat.print("kl block");
              }

            }
          }

          // contract bra and ket blocks
          const size_t nbra_ij = bra_ij.size();
          const size_t nket_ij = ket_ij.size();
          C_DGEMM('n', 't',
                  nbra_ij, nket_ij, blksize_int,
                  scale, bra_tile, blksize_int,
                  ket_tile, blksize_int,
                  0.0, T_result, tile_size_ket);
          if (debug >= DefaultPrintThresholds::allO2N2) {
            ExEnv::outn() << indent << "task " << me
                << ": nbra_ij = " << nbra_ij << " nket_ij = " << nket_ij << std::endl;
            ExEnv::outn() << indent << "task " << me
                << ": tile_size_bra = " << tile_size_bra << " tile_size_ket = " << tile_size_ket << std::endl;
            ExEnv::outn() << indent << "task " << me
                << ": T_result[0] = " << T_result[0] << std::endl;


            RefSCMatrix tmp_scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(tile_size_bra),
                                                                             new SCDimension(tile_size_ket));
            tmp_scmat.assign(T_result);
            tmp_scmat.print("ijkl");
          }

          // copy the result
          for(size_t ij=0; ij<nbra_ij; ++ij) {
            for(size_t kl=0; kl<nket_ij; ++kl) {

              const int i = bra_ij[ij].i0_;
              const int j = bra_ij[ij].i1_;
              const int k = ket_ij[kl].i0_;
              const int l = ket_ij[kl].i1_;

              //const double T_ijkl = T_result[ij * tile_size_ket + kl];
              braket->store_pair_subblock(i, j, 0,
                                          k, k+1, l, l+1,
                                          &(T_result[ij * tile_size_ket + kl]));
            }
          }

          if (ket_tile != 0) {
            for(size_t kl=0; kl<nket_ij; ++kl) {
              ket->release_pair_block(ket_ij[kl].i0_,
                                         ket_ij[kl].i1_,
                                         intsetidx_ket);
            }
            ket_tile = 0;
            ket_ij.resize(0);
          }

        } // ket tile loop

        if (bra_tile != 0) {
          const size_t nbra_ij = bra_ij.size();
          for(size_t ij=0; ij<nbra_ij; ++ij) {
            bra->release_pair_block(bra_ij[ij].i0_,
                                       bra_ij[ij].i1_,
                                       intsetidx_bra);
          }
          bra_tile = 0;
          bra_ij.resize(0);
        }

      } // bra tile loop
    } // loop over tasks with access

    if (bra->data_persistent()) bra->deactivate();
    if (bra != ket && ket->data_persistent()) ket->deactivate();
    if (braket->data_persistent()) braket->deactivate();

  }

  RefSCMatrix&
  operator<<(RefSCMatrix& dst,
             const Ref<DistArray4>& src) {

    assert(src->num_te_types() == 1);

    // is dst bra packed?
    bool bra_packed = false;
    const size_t nbra_sq = src->ni() * src->nj();
    const size_t nbra_tri = src->ni() * (src->ni() - 1) / 2;
    if (src->ni() != src->nj()) // no
      assert(dst.nrow() == nbra_sq);
    else { // maybe
      if (dst.nrow() == nbra_tri)
        bra_packed = true;
      else
        assert(dst.nrow() == nbra_sq);
    }

    // is dst ket packed?
    bool ket_packed = false;
    const size_t nket_sq = src->nx() * src->ny();
    const size_t nket_tri = src->nx() * (src->nx() - 1) / 2;
    if (src->nx() != src->ny()) // no
      assert(dst.ncol() == nket_sq);
    else { // maybe
      if (dst.ncol() == nket_tri)
        ket_packed = true;
      else
        assert(dst.ncol() == nket_sq);
    }

    RefSCVector row = dst.kit()->vector(dst.coldim());

    src->activate();
    SpinMOPairIter bra_iter(src->ni(), src->nj(), bra_packed);
    for(bra_iter.start(); bra_iter; bra_iter.next()) {
      const unsigned int b1 = bra_iter.i();
      const unsigned int b2 = bra_iter.j();
      const unsigned int b12 = bra_iter.ij();
      const double* b12_blk = src->retrieve_pair_block(b1, b2, 0);

      if (!ket_packed) {
        row.assign(b12_blk);
      }
      else {
        SpinMOPairIter ket_iter(src->nx(), src->ny(), ket_packed);
        for(ket_iter.start(); ket_iter; ket_iter.next()) {
          const unsigned int k1 = ket_iter.i();
          const unsigned int k2 = ket_iter.j();
          const unsigned int k12 = ket_iter.ij();
          row[k12] = b12_blk[k1 * src->ny() + k2];
        }
      }
      dst.assign_row(row, b12);

      src->release_pair_block(b1, b2, 0);

    }

    if (src->data_persistent()) src->deactivate();
  }

} // end of namespace sc

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
