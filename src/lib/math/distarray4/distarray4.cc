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

#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/misc/consumableresources.h>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>
#include <math/scmat/local.h>
#include <math/distarray4/distarray4.h>
#include <math/mmisc/pairiter.h>
#include <math/scmat/blas.h>
#include <util/misc/print.h>
#include <math/distarray4/distarray4_memgrp.h>
#include <math/distarray4/distarray4_node0file.h>
#ifdef HAVE_MPI
#include <math/distarray4/distarray4_mpiiofile.h>
#endif

using namespace std;
using namespace sc;

namespace {

  // replace with standard tuple when we switch C++0x
  template <typename T0, typename T1, typename T2>
  struct triple {
    public:
      triple() {}
      triple(const T0& i0,
             const T1& i1,
             const T2& i2) : i0_(i0), i1_(i1), i2_(i2) {}
      triple(const triple& other) : i0_(other.i0_), i1_(other.i1_), i2_(other.i2_) {}
      triple& operator=(const triple& other) {
        i0_ = other.i0_;
        i1_ = other.i1_;
        i2_ = other.i2_;
        return *this;
      }
      T0 i0_;
      T1 i1_;
      T2 i2_;
  };

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

  twa_map.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (has_access(proc)) {
      twa_map[proc] = count;
      count++;
    }
    else
      twa_map[proc] = -1;

  return count;
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
  Ref<DistArray4> permute23(const Ref<DistArray4>& src) {

    const size_t nt = src->num_te_types();
    const size_t ni = src->ni();
    const size_t nj = src->nj();
    const size_t nx = src->nx();
    const size_t ny = src->ny();
    const size_t njy = nj * ny;
    DistArray4Dimensions result_dims(nt, ni, nx, nj, ny, DistArray4Storage_XY);
    Ref<DistArray4> result = src->clone(result_dims);
    src->activate();
    result->activate();

    // determine the size and number of x tiles
    const size_t jy_blksize = njy * sizeof(double);
    const size_t available_memory = ConsumableResources::get_default_instance()->memory();
    size_t tilesize = (available_memory + jy_blksize - 1) / jy_blksize;
    if (tilesize > nx)  tilesize = nx;
    const size_t ntiles = (nx + tilesize - 1) / tilesize;
    tilesize = (nx + ntiles - 1) / ntiles;  // recompute the tile size for a given number of tiles
                                            // to make the tiles more even

    // allocate memory for the work buffers
    double* xjy_buf = allocate<double>(tilesize * njy);
    double* xy_buf = allocate<double>(tilesize * ny);

    // determine how many workers we have
    std::vector<int> worker_id;
    const int nworkers = src->tasks_with_access(worker_id);
    const int me = src->msg()->me();

    // sort!
    int task_id = 0;
    for(int i=0; i<ni; ++i) {
      for(int t=0; t<nt; ++t) {

        size_t xoffset = 0;
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

    deallocate(xjy_buf);
    deallocate(xy_buf);
    if(src->data_persistent()) src->deactivate();
    if(result->data_persistent()) result->deactivate();

    return result;
  }



  Ref<DistArray4> permute34(const Ref<DistArray4>& src)
  {
    const int ni = src->ni();
    const int nj = src->nj();
    const int nx = src->nx();
    const int ny = src->ny();
    const unsigned int ntypes = src->num_te_types();
    DistArray4Dimensions result_dims(ntypes, ni, nj, ny, nx, DistArray4Storage_XY);
    Ref<DistArray4> result = src->clone(result_dims);
    src->activate();
    result->activate();

    double* tmp_blk = allocate<double>(nx * ny);

    for (int t = 0; t < ntypes; ++t)
    {
      for (unsigned int b1 = 0; b1 < ni; ++b1)
      {
        for (unsigned int b2 = 0; b2 < nj; ++b2)
        {
          const double * b1b2_blk = src->retrieve_pair_block(b1, b2, t);
          for (unsigned int k1 = 0; k1 < nx; ++k1)
          {
            for (unsigned int k2 = 0; k2 < ny; ++k2)
            {
              tmp_blk[k1 * ny + k2] = b1b2_blk[k2 * nx + k1];
            }
          }
          result->store_pair_block(b1, b2, t, tmp_blk);
          src->release_pair_block(b1, b2, t);
        }
      }
    }

    deallocate(tmp_blk);
    if(src->data_persistent()) src->deactivate();
    if(result->data_persistent()) result->deactivate();
    return result;
  }


  Ref<DistArray4> permute12(const Ref<DistArray4>& src)
 {
    const int ni = src->ni();
    const int nj = src->nj();
    const int nx = src->nx();
    const int ny = src->ny();
    const unsigned int ntypes = src->num_te_types();
    DistArray4Dimensions result_dims(ntypes, nj, ni, nx, ny, DistArray4Storage_XY);
    Ref<DistArray4> result = src->clone(result_dims);
    src->activate();
    result->activate();

    double* tmp_blk = allocate<double>(nx * ny);

    for (int t = 0; t < ntypes; ++t)
    {
      for (unsigned int b1 = 0; b1 < ni; ++b1)
      {
        for (unsigned int b2 = 0; b2 < nj; ++b2)
        {
          const double * b1b2_blk = src->retrieve_pair_block(b1, b2, t, tmp_blk);
          result->store_pair_block(b2, b1, t, tmp_blk);
          src->release_pair_block(b1, b2, t);

        }
      }
    }

    deallocate(tmp_blk);
    if(src->data_persistent()) src->deactivate();
    if(result->data_persistent()) result->deactivate();
    return result;
  }


  void axpy(const Ref<DistArray4>& X,
            double a,
            const Ref<DistArray4>& Y,
            double scale) {

    MPQC_ASSERT( X->num_te_types() == Y->num_te_types() );
    MPQC_ASSERT( X->ni() == Y->ni() );
    MPQC_ASSERT( X->nj() == Y->nj() );
    MPQC_ASSERT( X->nx() == Y->nx() );
    MPQC_ASSERT( X->ny() == Y->ny() );
    const int nt = Y->num_te_types();
    const int ni = Y->ni();
    const int nj = Y->nj();
    const int nx = Y->nx();
    const int ny = Y->ny();
    const blasint nxy = nx * ny;
    const size_t bufsize = nxy * sizeof(double);

    X->activate();
    Y->activate();

    // maximum tile size is determined by the available memory
    double* y_buf = allocate<double>(nxy);

    // determine how many workers we have
    std::vector<int> worker_id;
    const int nworkers = Y->tasks_with_access(worker_id);
    const int me = Y->msg()->me();

    size_t task_id = 0;
    for(int i=0; i<ni; ++i) {
      for(int j=0; j<nj; ++j) {
        for(int t=0; t<nt; ++t, ++task_id) {

          // round-robin task allocation
          if(task_id%nworkers != worker_id[me])
            continue;

          const double* x_buf = X->retrieve_pair_block(i, j, t);
          Y->retrieve_pair_block(i, j, t, y_buf);

          const blasint one = 1;
          F77_DAXPY(&nxy, &a, x_buf, &one, y_buf, &one);
          if (scale != 1.0) F77_DSCAL(&nxy, &scale, y_buf, &one);

          Y->store_pair_block(i, j, t, y_buf);

          X->release_pair_block(i, j, t);
          Y->release_pair_block(i, j, t);
        }

      }
    }

    deallocate(y_buf);

    if (X->data_persistent()) X->deactivate();
    if (Y->data_persistent()) Y->deactivate();
  }




  void contract34(Ref<DistArray4>& braket,
                  double scale,
                  const Ref<DistArray4>& bra,
                  unsigned int intsetidx_bra,
                  const Ref<DistArray4>& ket,
                  unsigned int intsetidx_ket,
                  int debug) {
    if(braket == 0)
    {
      DistArray4Dimensions  braket_dims(1, bra->ni(), bra->nj(),
                                        ket->ni(), ket->nj());
      braket = bra->clone(braket_dims);
    }

    bra->activate();
    ket->activate();
    braket->activate();

    const unsigned int nb1 = bra->ni();
    const unsigned int nb2 = bra->nj();
    const unsigned int nk1 = ket->ni();
    const unsigned int nk2 = ket->nj();
    const unsigned int n1 =  bra->nx();
    const unsigned int n2 =  bra->ny();
    MPQC_ASSERT(n1 == ket->nx());
    MPQC_ASSERT(n2 == ket->ny());
    MPQC_ASSERT(braket->ni() == nb1);
    MPQC_ASSERT(braket->nj() == nb2);
    MPQC_ASSERT(braket->nx() == nk1);
    MPQC_ASSERT(braket->ny() == nk2);

    // Using spinorbital iterators means I don't take into account perm symmetry
    // More efficient algorithm will require generic code
    typedef sc::fastpairiter::MOPairIter<sc::fastpairiter::ASymm> MOPairIter;
    MOPairIter iterbra(nb1, nb2);
    MOPairIter iterket(nk1, nk2);
    MOPairIter iterint( n1,  n2);
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
    const size_t nproc = bra->msg()->n();
    size_t try_tile_size_bra = (nij_bra + nproc - 1) / nproc;
    size_t try_tile_size_ket = (nij_ket + nproc - 1) / nproc;
    try_tile_size_bra = std::min(try_tile_size_bra, max_tile_size);
    try_tile_size_ket = std::min(try_tile_size_ket, max_tile_size);
    const size_t ntiles_bra = (nij_bra + try_tile_size_bra - 1) / try_tile_size_bra;
    const size_t ntiles_ket = (nij_ket + try_tile_size_ket - 1) / try_tile_size_ket;
    const size_t tile_size_bra = (nij_bra + ntiles_bra - 1) / ntiles_bra;
    const size_t tile_size_ket = (nij_ket + ntiles_ket - 1) / ntiles_ket;

    // scratch buffers to hold tiles and the result of the contraction
    double* T_bra = allocate<double>(tile_size_bra * blksize_int);
    double* T_ket = allocate<double>(tile_size_ket * blksize_int);
    double* T_result = allocate<double>(tile_size_bra * tile_size_ket);

    // split work over tasks which have access to integrals
    // WARNING: assuming same accessibility for both bra and ket transforms
    std::vector<int> proc_with_ints;
    const int nproc_with_ints = bra->tasks_with_access(proc_with_ints);
    const int me = bra->msg()->me();

    if (bra->has_access(me)) {

      size_t task_count = 0;

      // these will keep track of each tile's sets of i and j values
      typedef triple<unsigned int, unsigned int, unsigned int> uint3;
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

    deallocate(T_bra);
    deallocate(T_ket);
    deallocate(T_result);

    if (bra->data_persistent()) bra->deactivate();
    if (bra != ket && ket->data_persistent()) ket->deactivate();
    if (braket->data_persistent()) braket->deactivate();

    bra->msg()->sync();
  }



  void contract34_DA4_RefMat(Ref<DistArray4>& braket,
                                  double scale,
                                  const Ref<DistArray4>& bra,
                                  unsigned int intsetidx_bra,
                                  const RefSCMatrix& ket,
                                  const int MatBra1Dim, const int MatBra2Dim)
  {
    if(braket == 0)
    {
      DistArray4Dimensions  braket_dims(1, bra->ni(), bra->nj(),
                                           MatBra1Dim, MatBra2Dim);
      braket = bra->clone(braket_dims);
    }
    bra->activate();
    braket->activate();

    const unsigned int nb1 = bra->ni();
    const unsigned int nb2 = bra->nj();
    const unsigned int nk1 = MatBra1Dim;
    const unsigned int nk2 = MatBra2Dim;
    const unsigned int n1 =  bra->nx();
    const unsigned int n2 =  bra->ny();
    MPQC_ASSERT(braket->ni() == nb1);
    MPQC_ASSERT(braket->nj() == nb2);
    MPQC_ASSERT(braket->nx() == MatBra1Dim);
    MPQC_ASSERT(braket->ny() == MatBra2Dim);
//    abort();

    // Using spinorbital iterators means I don't take into account perm symmetry
    // More efficient algorithm will require generic code
    typedef sc::fastpairiter::MOPairIter<sc::fastpairiter::ASymm> MOPairIter;
    MOPairIter iterbra(nb1, nb2);
    MOPairIter iterket(nk1, nk2);
    MOPairIter iterint(n1,  n2);
    // size of one block of <space1_bra space2_bra|
    const unsigned int nbra = iterbra.nij();
    // size of one block of <space1_ket space2_ket|
    const unsigned int nket = iterket.nij();

    // data will be accessed in tiles
    // determine tiling here
    //
    // size of each integral block |int1 int2)
    const unsigned int blksize_int = n1 * n2;


    // maximum tile size is determined by the available memory
    const size_t memory_available = ConsumableResources::get_default_instance()->memory();
    const size_t max_tile_size = (memory_available - nket *blksize_int*sizeof(double)) / (blksize_int * sizeof(double));
    if (max_tile_size == 0) {
      throw AlgorithmException("not enough memory for a single tile, increase memory", __FILE__, __LINE__);
    }
    // try tiling nij_bra into nproc tiles
    const int nproc = bra->msg()->n();
    size_t try_tile_size_bra = (nbra + nproc - 1) / nproc;
    try_tile_size_bra = std::min(try_tile_size_bra, max_tile_size);
    const size_t ntiles_bra = (nbra + try_tile_size_bra - 1) / try_tile_size_bra;
    const size_t tile_size_bra = (nbra + ntiles_bra - 1) / ntiles_bra;

    // scratch buffers to hold tiles and the result of the contraction
    double* T_bra = allocate<double>(tile_size_bra * blksize_int);
    double* T_ket = allocate<double>(nket * blksize_int);
    ket.convert(T_ket);
    double* T_result = allocate<double>(tile_size_bra * nket);


    // split work over tasks which have access to integrals
    // WARNING: assuming same accessibility for both bra and ket transforms
    std::vector<int> proc_with_ints;
    const int nproc_with_ints = bra->tasks_with_access(proc_with_ints);
    const int me = bra->msg()->me();

    if (bra->has_access(me)) {

      size_t task_count = 0; // task counts starts

      typedef triple<unsigned int, unsigned int, unsigned int> uint3; // keep track of each tile's sets of i and j values
      std::vector<uint3> bra_ij;
      iterbra.start();
      size_t tbra_offset = 0; // 'tbra' loops over bra tiles for this set; 'tbra_offset' tracks the block pos after each tile
      for (size_t tbra = 0; tbra < ntiles_bra; ++tbra, tbra_offset += tile_size_bra, ++task_count) {
        double* bra_tile = 0;   // zero indicates it needs to be loaded

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

              // where to read the integrals? if not antisymmetrizing read directly to T_bra, else read to scratch
              // buffer, then antisymmetrize to T_bra
              double* blk_ptr_read = bra_tile + i*blksize_int;
              Timer tim_intsretrieve("MO ints retrieve");
              const double* blk_cptr = bra->retrieve_pair_block(ijt.i0_,
                                                                   ijt.i1_,
                                                                   intsetidx_bra,
                                                                   blk_ptr_read);
              tim_intsretrieve.exit();
            }
          }


          // contract bra and ket blocks
          const size_t nbra_ij = bra_ij.size();
          C_DGEMM('n', 't',
                  nbra_ij, nket, blksize_int,
                  scale, bra_tile, blksize_int,
                  T_ket, blksize_int,
                  0.0, T_result, nket);


          // copy the result
          for(size_t ij=0; ij<nbra_ij; ++ij) {
              const int i = bra_ij[ij].i0_;
              const int j = bra_ij[ij].i1_;
              braket->store_pair_block(i, j, 0, &(T_result[ij * nket]));
            }

        if (bra_tile != 0) {
          const size_t nbra_ij = bra_ij.size();
          for(size_t ij=0; ij<nbra_ij; ++ij)
          {
            bra->release_pair_block(bra_ij[ij].i0_,
                                       bra_ij[ij].i1_,
                                       intsetidx_bra);
          }
          bra_tile = 0;
          bra_ij.resize(0);
        }

      } // bra tile loop
    } // loop over tasks with access

    deallocate(T_bra);
    deallocate(T_ket);
    deallocate(T_result);

    if (bra->data_persistent()) bra->deactivate();
    if (braket->data_persistent()) braket->deactivate();
  }


  void contract_DA4_RefMat_k2b2_34(Ref<DistArray4>& braket,
                         double scale,
                         const Ref<DistArray4>& bra,
                         unsigned int intsetidx_bra,
                         const RefSCMatrix& ket,
                         const int MatBra1Dim, const int MatBra2Dim)
  {
    if(braket == 0)
     {
       DistArray4Dimensions  braket_dims(1, bra->ni(), MatBra1Dim, bra->nx(), MatBra2Dim);
       braket = bra->clone(braket_dims);
     }
    Ref<DistArray4> DA_Aixy = permute23(bra);
    Ref<DistArray4> DA_M;
    contract34_DA4_RefMat(DA_M, scale, DA_Aixy,intsetidx_bra, ket, MatBra1Dim, MatBra2Dim);
    braket = permute23(DA_M);
    return;
  }

  void contract_DA4_RefMat_k1b2_34(Ref<DistArray4>& braket,
                         double scale,
                         const Ref<DistArray4>& bra,
                         unsigned int intsetidx_bra,
                         const RefSCMatrix& ket,
                         const int MatBra1Dim, const int MatBra2Dim)
  {
    if(braket == 0)
     {
       DistArray4Dimensions  braket_dims(1, bra->ni(), MatBra1Dim, bra->ny(), MatBra2Dim);
       braket = bra->clone(braket_dims);
     }

    Ref<DistArray4> DA_Aixy = permute23(permute34(bra));
    Ref<DistArray4> DA_M;
    contract34_DA4_RefMat(DA_M, scale, DA_Aixy,intsetidx_bra, ket, MatBra1Dim, MatBra2Dim);
    braket = permute23(DA_M);
    return;
  }


  namespace {

    enum Index34 {
      Index3, Index4
    };

    /// computes (ij|xy) . T_x^z  or same for index 4
    /// note that the result is not symmetric with respect to permutation of particles 1 and 2
    template <Index34 ContrIndex>
    void
    _contract_3_or_4(const Ref<DistArray4>& src, const RefSCMatrix& tform, Ref<DistArray4>& dest)
    {
      if (dest == 0) {
        DistArray4Dimensions dest_dims(src->num_te_types(), src->ni(), src->nj(),
                                       ((ContrIndex == Index3) ? tform.ncol() : src->nx()),
                                       ((ContrIndex == Index4) ? tform.ncol() : src->ny())
                                       );
        dest = src->clone(dest_dims);
      }

      MPQC_ASSERT(src->num_te_types() == dest->num_te_types());
      MPQC_ASSERT(src->ni() == dest->ni());
      MPQC_ASSERT(src->nj() == dest->nj());
      MPQC_ASSERT(src->nx() == ((ContrIndex == Index3) ? tform.nrow() : dest->nx()));
      MPQC_ASSERT(src->ny() == ((ContrIndex == Index4) ? tform.nrow() : dest->ny()));
      MPQC_ASSERT(dest->nx() == ((ContrIndex == Index3) ? tform.ncol() : src->nx()));
      MPQC_ASSERT(dest->ny() == ((ContrIndex == Index4) ? tform.ncol() : src->ny()));

      // copy T to an array
      double* tform_buf = allocate<double>(tform.nrow() * tform.ncol());
      tform.convert(tform_buf);

      // allocate space for the result
      double* dest_buf = allocate<double>(dest->nx() * dest->ny());

      // determine how many workers we have
      std::vector<int> worker_id;
      const int nworkers = dest->tasks_with_access(worker_id);
      const int me = dest->msg()->me();

      src->activate();
      dest->activate();

      const unsigned int ni = src->ni();
      const unsigned int nj = src->nj();
      const unsigned int nx = src->nx();
      const unsigned int ny = src->ny();
      const unsigned int nX = dest->nx();
      const unsigned int nY = dest->ny();
      size_t task_id = 0;
      for(int i=0; i<ni; ++i) {
        for(int j=0; j<nj; ++j) {
          for(int t=0; t<dest->num_te_types(); ++t, ++task_id) {

            // round-robin task allocation
            if(task_id%nworkers != worker_id[me])
              continue;

            const double* src_buf = src->retrieve_pair_block(i, j, t);

            if (ContrIndex == Index4)
              // src_buf * tform_buf = dest_buf
              C_DGEMM('n','n', nx, nY, ny, 1.0, src_buf, ny, tform_buf, nY, 0.0, dest_buf, nY);
            else  // ContrIndex == Index3
              // tform_buf^t * src_buf  = dest_buf
              C_DGEMM('t','n', nX, nY, nx, 1.0, tform_buf, nX, src_buf, ny, 0.0, dest_buf, nY);

            dest->store_pair_block(i, j, t, dest_buf);
            src->release_pair_block(i, j, t);
          }
        }
      }

      if (src->data_persistent()) src->deactivate();
      if (dest->data_persistent()) dest->deactivate();

      deallocate(dest_buf);
      deallocate(tform_buf);
    }

    template <bool Antisymmetrize>
    void
    _symmetrize(const Ref<DistArray4>& A)
    {
      MPQC_ASSERT(A->ni() == A->nj());
      MPQC_ASSERT(A->nx() == A->ny());

      A->activate();

      const unsigned int nbra = A->ni();
      const unsigned int nket = A->nx();
      const unsigned int ntypes = A->num_te_types();
      double* tmp_blk = allocate<double>(nket * nket);

      for (int t = 0; t < ntypes; ++t) {
        for (unsigned int b1 = 0; b1 < nbra; ++b1) {
          for (unsigned int b2 = 0; b2 <= b1; ++b2) {

            const double* b12_blk = A->retrieve_pair_block(b1, b2, t);
            const double* b21_blk = A->retrieve_pair_block(b2, b1, t);

            size_t k12 = 0;
            for (unsigned int k1 = 0; k1 < nket; ++k1) {
              size_t k21 = k1;
              for (unsigned int k2 = 0; k2 < nket; ++k2, ++k12, k21 += nket) {

                if (Antisymmetrize)
                  tmp_blk[k12] = 0.5 * (b12_blk[k12] + b21_blk[k21] - b21_blk[k12] - b12_blk[k21]);
                else
                  tmp_blk[k12] = 0.5 * (b12_blk[k12] + b21_blk[k21]);

              }
            }

            A->release_pair_block(b1, b2, t);
            A->release_pair_block(b2, b1, t);

            A->store_pair_block(b1, b2, t, tmp_blk);
            if (b1 != b2) {
              for (unsigned int k1 = 0; k1 < nket; ++k1) {
                for (unsigned int k2 = 0; k2 <= k1; ++k2) {
                  const size_t k12 = k1 * nket + k2;
                  const size_t k21 = k2 * nket + k1;
                  const double tmp = tmp_blk[k12];
                  tmp_blk[k12] = tmp_blk[k21];
                  tmp_blk[k21] = tmp;
                }
              }
              A->store_pair_block(b2, b1, t, tmp_blk);
            }
          }
        }
      }

      deallocate(tmp_blk);
      if (A->data_persistent()) A->deactivate();
    }

  } // end of sc::anonymous namespace

  void contract3(const Ref<DistArray4>& ijxy, const RefSCMatrix& T, Ref<DistArray4>& ijzy){
    _contract_3_or_4<Index3>(ijxy, T, ijzy);
  }
  void contract4(const Ref<DistArray4>& ijxy, const RefSCMatrix& T, Ref<DistArray4>& ijxz){
    _contract_3_or_4<Index4>(ijxy, T, ijxz);
  }

  void antisymmetrize(const Ref<DistArray4>& A) {
    _symmetrize<true>(A);
  }
  void symmetrize(const Ref<DistArray4>& A) {
    _symmetrize<false>(A);
  }

  Ref<DistArray4> extract(const Ref<DistArray4>& A,
                          unsigned int te_type,
                          double scale) {
    MPQC_ASSERT(te_type < A->num_te_types());

    DistArray4Dimensions dims(1, A->ni(), A->nj(), A->nx(), A->ny());
    Ref<DistArray4> result = A->clone(dims);

    A->activate();
    result->activate();

    // determine how many workers we have
    std::vector<int> worker_id;
    const int nworkers = A->tasks_with_access(worker_id);
    const int me = A->msg()->me();

    const blasint nxy = A->nx() * A->ny();

    for (unsigned int i = 0, task_id = 0; i < A->ni(); ++i) {
      for (unsigned int j = 0; j < A->nj(); ++j, ++task_id) {

        // round-robin task allocation
        if(task_id%nworkers != worker_id[me])
          continue;

        const double* ij_blk = A->retrieve_pair_block(i, j, te_type);
        if (scale != 1.0) {
          const blasint one = 1;
          F77_DSCAL(&nxy, &scale, const_cast<double*>(ij_blk), &one);
        }
        result->store_pair_block(i, j, 0, ij_blk);
        A->release_pair_block(i, j, te_type);
      }
    }

    if (A->data_persistent()) A->deactivate();
    if (result->data_persistent()) result->deactivate();
    return result;
  }


  RefSCMatrix &
  copy_to_RefSCMat(RefSCMatrix& dst,
                          const Ref<DistArray4>& src, const int tensor_index)
  {
    MPQC_ASSERT(src->has_access(src->msg()->me()));

    // is dst bra packed?
    bool bra_packed = false;
    const size_t nbra_sq = src->ni() * src->nj();
    const size_t nbra_tri = src->ni() * (src->ni() - 1) / 2;
    if (src->ni() != src->nj()) // no
      MPQC_ASSERT(dst.nrow() == nbra_sq);
    else { // maybe
      if (dst.nrow() == nbra_tri)
        bra_packed = true;
      else
        MPQC_ASSERT(dst.nrow() == nbra_sq);
    }

    // is dst ket packed?
    bool ket_packed = false;
    const size_t nket_sq = src->nx() * src->ny();
    const size_t nket_tri = src->nx() * (src->nx() - 1) / 2;
    if (src->nx() != src->ny()) // no
      MPQC_ASSERT(dst.ncol() == nket_sq);
    else { // maybe
      if (dst.ncol() == nket_tri)
        ket_packed = true;
      else
        MPQC_ASSERT(dst.ncol() == nket_sq);
    }

    RefSCVector row = dst.kit()->vector(dst.coldim());

    src->activate();
    SpinMOPairIter bra_iter(src->ni(), src->nj(), bra_packed);
    for(bra_iter.start(); bra_iter; bra_iter.next()) {
      const unsigned int b1 = bra_iter.i();
      const unsigned int b2 = bra_iter.j();
      const unsigned int b12 = bra_iter.ij();
      const double* b12_blk = src->retrieve_pair_block(b1, b2, tensor_index);

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

      src->release_pair_block(b1, b2, tensor_index);
    }

    if (src->data_persistent()) src->deactivate();
    return dst;
  }

  RefSCMatrix
  copy_to_RefSCMat(const Ref<DistArray4>& src, int tensor_index)
  {
    RefSCDimension upp_scdim = new SCDimension(src->ni() * src->nj());
    RefSCDimension low_scdim = new SCDimension(src->nx() * src->ny());
    Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
    RefSCMatrix dst = local_kit->matrix(upp_scdim, low_scdim);
    dst.assign(0.0);
    copy_to_RefSCMat(dst, src, tensor_index);
    return dst;
  }

 RefSCMatrix &
  operator<<(RefSCMatrix& dst, const Ref<DistArray4>& src)
  {
    MPQC_ASSERT(src->num_te_types() == 1);
    copy_to_RefSCMat(dst, src, 0);
    return dst;
  }



 Ref<DistArray4> make_distarray4(int num_te_types, int ni, int nj, int nx, int ny,
                                 DistArray4Storage storage) {

   const int nproc = MessageGrp::get_default_messagegrp()->n();
   const size_t nij_local_max = ((size_t)ni*nj + nproc - 1)/nproc;
   const size_t blksize = num_te_types * nx * ny * sizeof(double);
   const size_t max_local_memory = nij_local_max * blksize;

   if (false && ConsumableResources::get_default_instance()->memory() > max_local_memory*3) {
     return new DistArray4_MemoryGrp(MemoryGrp::get_default_memorygrp(), num_te_types, ni, nj, nx, ny, blksize, storage);
   }
   else {
     const std::string basename = SCFormIO::fileext_to_filename(".moints");
     const std::string dir = ConsumableResources::get_default_instance()->disk_location();
     const std::string fname(tempnam(dir.c_str(), basename.c_str()));
//#if HAVE_MPIIO
//     return new DistArray4_MPIIOFile_Ind(fname.c_str(), num_te_types, ni, nj, nx, ny, storage);
//#else
     return new DistArray4_Node0File(fname.c_str(), num_te_types, ni, nj, nx, ny, storage);
//#endif
   }


 }



} // end of namespace sc

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
