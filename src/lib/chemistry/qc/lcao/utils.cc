//
// utils.cc
//
// Copyright (C) 2005 Edward Valeev
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
#include <util/misc/scexception.h>
#include <util/misc/consumableresources.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/qc/lcao/utils.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/wfn/orbitalspace.h>

using namespace sc;

void
sc::antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                   const Ref<OrbitalSpace>& bra,
                   const Ref<OrbitalSpace>& ket,
                   bool accumulate)
{
  SpatialMOPairIter_eq ij_iter(bra->rank());
  SpatialMOPairIter_eq kl_iter(ket->rank());
  const unsigned int brablock_size_ab = ij_iter.nij_ab();
  const unsigned int ketblock_size_ab = kl_iter.nij_ab();
  if (brablock_size_ab==0 || ketblock_size_ab==0)
    return;
  if (A.rowdim().n()%brablock_size_ab)
    throw ProgrammingError("sc::antisymmetrize() -- row dimension is not integer multiple of bra-space rank",__FILE__,__LINE__);
  if (A.coldim().n()%ketblock_size_ab)
    throw ProgrammingError("sc::antisymmetrize() -- col dimension is not integer multiple of ket-space rank",__FILE__,__LINE__);
  const unsigned int nbra_blocks = A.rowdim().n() / brablock_size_ab;
  const unsigned int nket_blocks = A.coldim().n() / ketblock_size_ab;
  const unsigned int brablock_size_aa = ij_iter.nij_aa();
  const unsigned int ketblock_size_aa = kl_iter.nij_aa();
  
  unsigned int bra_offset_ab = 0;
  unsigned int bra_offset_aa = 0;
  for(int brablock=0; brablock<nbra_blocks; brablock++, bra_offset_ab += brablock_size_ab, bra_offset_aa += brablock_size_aa) {
    for(ij_iter.start();int(ij_iter);ij_iter.next()) {
      
      const int ij_aa = ij_iter.ij_aa();
      if (ij_aa == -1)
        continue;
      const int ij_ab = ij_iter.ij_ab();
      
      unsigned int ket_offset_ab = 0;
      unsigned int ket_offset_aa = 0;
      for(int ketblock=0; ketblock<nket_blocks; ketblock++, ket_offset_ab += ketblock_size_ab, ket_offset_aa += ketblock_size_aa) {
        for(kl_iter.start();int(kl_iter);kl_iter.next()) {
          
          const int kl_aa = kl_iter.ij_aa();
          if (kl_aa == -1)
            continue;
          const int kl_ab = kl_iter.ij_ab();
          const int lk_ab = kl_iter.ij_ba();
          
          double Aanti_ijkl = A.get_element(ij_ab+bra_offset_ab,kl_ab+ket_offset_ab) - A.get_element(ij_ab+bra_offset_ab,lk_ab+ket_offset_ab);
          if (accumulate)
            Aanti.accumulate_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
          else
            Aanti.set_element(ij_aa+bra_offset_aa,kl_aa+ket_offset_aa,Aanti_ijkl);
        }
      }
    }
  }
}

std::vector<double>
sc::convert(const RefDiagSCMatrix& A)
{
  const int n = A.dim().n();
  std::vector<double> result;
  for(int i=0; i<n; i++)
    result.push_back(A.get_element(i));
  return result;
}

void
sc::print_f77_mat(const std::string& comment,
                  const double* A,
                  unsigned int nrow,
                  unsigned int ncol,
                  bool transpose)
{
  RefSCDimension rowdim = new SCDimension(nrow);
  RefSCDimension coldim = new SCDimension(ncol);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix amat = localkit->matrix(rowdim,coldim);
  amat.assign(A);
  if (transpose) amat = amat.t();
  amat.print(comment.c_str());
}

/// Returns the lower triangle of the matrix B (which should be symmetric)
sc::RefSymmSCMatrix
sc::to_lower_triangle(const RefSCMatrix& B) {
  RefSymmSCMatrix Bs = B.kit()->symmmatrix(B.rowdim());
  int n = B.nrow();
  double* b = new double[n*n];
  B.convert(b);
  const double* b_ptr = b;
  for(int i=0; i<n; i++, b_ptr += i)
    for(int j=i; j<n; j++, b_ptr++)
      Bs.set_element(i,j,*b_ptr);
  delete[] b;
  return Bs;
}

void
sc::map(const Ref<DistArray4>& src,
        const Ref<OrbitalSpace>& isrc,
        const Ref<OrbitalSpace>& jsrc,
        const Ref<OrbitalSpace>& xsrc,
        const Ref<OrbitalSpace>& ysrc,
        Ref<DistArray4>& dest,
        const Ref<OrbitalSpace>& idest,
        const Ref<OrbitalSpace>& jdest,
        const Ref<OrbitalSpace>& xdest,
        const Ref<OrbitalSpace>& ydest) {

  if (dest == 0) { // allocate dest, if needed, by cloning src
    DistArray4Dimensions dest_dims(src->num_te_types(), idest->rank(), jdest->rank(), xdest->rank(), ydest->rank());
    dest = src->clone(dest_dims);
  }

  MPQC_ASSERT( src->num_te_types() == dest->num_te_types() );
  MPQC_ASSERT( src->ni() == isrc->rank() );
  MPQC_ASSERT( src->nj() == jsrc->rank() );
  MPQC_ASSERT( src->nx() == xsrc->rank() );
  MPQC_ASSERT( src->ny() == ysrc->rank() );
  MPQC_ASSERT( dest->ni() == idest->rank() );
  MPQC_ASSERT( dest->nj() == jdest->rank() );
  MPQC_ASSERT( dest->nx() == xdest->rank() );
  MPQC_ASSERT( dest->ny() == ydest->rank() );

  MOIndexMap imap(*isrc << *idest);
  MOIndexMap jmap(*jsrc << *jdest);
  MOIndexMap xmap(*xsrc << *xdest);
  MOIndexMap ymap(*ysrc << *ydest);

  src->activate();
  dest->activate();

  double* dest_buf = allocate<double>(dest->nx() * dest->ny());

  // determine how many workers we have
  std::vector<int> worker_id;
  const int nworkers = dest->tasks_with_access(worker_id);
  const int me = dest->msg()->me();

  const unsigned int nx = xdest->rank();
  const unsigned int ny = ydest->rank();
  const unsigned int nX = xsrc->rank();
  const unsigned int nY = ysrc->rank();
  size_t task_id = 0;
  for(int i=0; i<idest->rank(); ++i) {
    for(int j=0; j<jdest->rank(); ++j) {
      for(int t=0; t<dest->num_te_types(); ++t, ++task_id) {

        // round-robin task allocation
        if(task_id%nworkers != worker_id[me])
          continue;

        const unsigned int I = imap[i];
        const unsigned int J = jmap[j];
        const double* src_buf = src->retrieve_pair_block(I, J, t);

        size_t xy = 0;
        for(unsigned int x=0; x<nx; ++x) {
          const unsigned int X = xmap[x];
          for(unsigned int y=0; y<ny; ++y, ++xy) {
            const unsigned int Y = ymap[y];
            const size_t XY = X * nY + Y;

            dest_buf[xy] = src_buf[XY];
          }
        }
        dest->store_pair_block(i, j, t, dest_buf);
        src->release_pair_block(I, J, t);

      }
    }
  }

  deallocate(dest_buf);

  if (src->data_persistent()) src->deactivate();
  if (dest->data_persistent()) dest->deactivate();
}

///////
Ref<sc::SpatialMOPairIter>
MOPairIterFactory::mopairiter(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2)
{
  if (space1 == space2)
    return new SpatialMOPairIter_eq(space1->rank());
  else
    return new SpatialMOPairIter_neq(space1->rank(), space2->rank());
}

RefSCDimension
MOPairIterFactory::scdim_aa(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2)
{
  if (space1 != space2)
    return scdim_ab(space1,space2);
  else {
    const int n = space1->rank();
    const int npair_aa = n*(n-1)/2;
    std::string name = "Alpha-alpha pair (" + space1->name() + "," + space1->name() + ")";
    return new SCDimension(npair_aa,name.c_str());
  }
}

RefSCDimension
MOPairIterFactory::scdim_ab(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2)
{
  Ref<MOPairIter> piter = mopairiter(space1,space2);
  int npair_ab = space1->rank() * space2->rank();
  std::string name = "Alpha-beta pair (" + space1->name() + "," + space2->name() + ")";
  return new SCDimension(npair_ab,name.c_str());
}
