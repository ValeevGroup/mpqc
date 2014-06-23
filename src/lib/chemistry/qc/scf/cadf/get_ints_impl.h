//
// get_ints_impl.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Apr 18, 2014
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

#ifndef _chemistry_qc_scf_get_ints_impl_h
#define _chemistry_qc_scf_get_ints_impl_h

#include "cadfclhf.h" // for the benefit of the parser

namespace sc {

template <typename ShellRange>
CADFCLHF::TwoCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<TwoCenterIntContainer>(iblk.nbf, jblk.nbf);
  int block_offset_i = 0;
  for(auto ish : shell_range(iblk)) {
    int block_offset_j = 0;
    for(auto jsh : shell_range(jblk)) {
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type, ish.basis, jsh.basis);
      rv->block(block_offset_i, block_offset_j, ish.nbf, jsh.nbf) = *ints_ptr;
      block_offset_j += jsh.nbf;
    }
    block_offset_i += ish.nbf;
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::TwoCenterIntContainerPtr
CADFCLHF::ints_to_eigen_threaded(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    std::vector<Ref<TwoBodyTwoCenterInt>>& ints_for_thread,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<TwoCenterIntContainer>(iblk.nbf, jblk.nbf);
  // non-contiguous form not implemented yet
  assert(iblk.is_contiguous());
  assert(jblk.is_contiguous());
  do_threaded(nthread_, [&](int ithr){
    ShellData ish, jsh;
    for(auto&& pair : threaded_shell_block_pair_range(iblk, jblk, ithr, nthread_)){
      boost::tie(ish, jsh) = pair;
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints_for_thread[ithr], int_type, ish.basis, jsh.basis);
      rv->block(
          ish.bfoff - iblk.bfoff, jsh.bfoff - jblk.bfoff,
          ish.nbf, jsh.nbf
      ) = *ints_ptr;
    }
  });
  return rv;
}

template <typename ShellRange>
Eigen::Map<CADFCLHF::TwoCenterIntContainer>
CADFCLHF::ints_to_eigen_map_threaded(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    std::vector<Ref<TwoBodyTwoCenterInt>>& ints_for_thread,
    TwoBodyOper::type int_type,
    double* buffer
){
  auto rv = Eigen::Map<TwoCenterIntContainer>(buffer, iblk.nbf, jblk.nbf);
  // non-contiguous form not implemented yet
  assert(iblk.is_contiguous());
  assert(jblk.is_contiguous());
  do_threaded(nthread_, [&](int ithr){
    ShellData ish, jsh;
    for(auto&& pair : threaded_shell_block_pair_range(iblk, jblk, ithr, nthread_)){
      boost::tie(ish, jsh) = pair;
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints_for_thread[ithr], int_type, ish.basis, jsh.basis);
      rv.block(
          ish.bfoff - iblk.bfoff, jsh.bfoff - jblk.bfoff,
          ish.nbf, jsh.nbf
      ) = *ints_ptr;
    }
  });
  return rv;
}


template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>((const int)iblk.nbf, (const int)jsh.nbf);
  for(auto ish : shell_range(iblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type);
    rv->block(ish.bfoff - iblk.bfoff, 0, ish.nbf, jsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellData& ish, const ShellData& jsh,
    const ShellBlockData<ShellRange>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      ish.nbf * jsh.nbf,
      Xblk.nbf
  );
  for(auto Xsh : shell_range(Xblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type, ish.basis, jsh.basis, Xblk.basis);
    rv->middleCols(Xsh.bfoff - Xblk.bfoff, Xsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange1, typename ShellRange2>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange1>& iblk,
    const ShellData& jsh,
    const ShellBlockData<ShellRange2>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      iblk.nbf * jsh.nbf,
      Xblk.nbf
  );
  for(auto ish: shell_range(iblk)) {
    for(auto Xsh : shell_range(Xblk)) {
      const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type, ish.basis, jsh.basis, Xsh.basis);
      rv->block(
          (ish.bfoff - iblk.bfoff) * jsh.nbf, Xsh.bfoff - Xblk.bfoff,
          ish.nbf*jsh.nbf, Xsh.nbf
      ) = *ints_ptr;
    }
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      iblk.nbf * jsh.nbf, Xsh.nbf
  );
  int block_offset = 0;
  for(auto ish: shell_range(iblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type, ish.basis, jsh.basis, Xsh.basis);
    rv->block(
        block_offset * jsh.nbf, 0,
        ish.nbf*jsh.nbf, Xsh.nbf
    ) = *ints_ptr;
    block_offset += ish.nbf;
  }
  return rv;
}

template <typename ShellRange>
Eigen::Map<CADFCLHF::ThreeCenterIntContainer>
CADFCLHF::ints_to_eigen_map(
    const ShellData& ish,
    const ShellData& jsh,
    const ShellBlockData<ShellRange>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type,
    double* __restrict__ buffer
){
  //
  Eigen::Map<ThreeCenterIntContainer> rv(
      buffer, ish.nbf * jsh.nbf, Xblk.nbf
  );
  int block_offset = 0;
  typedef Eigen::Map<RowMatrix, Eigen::Default, Eigen::OuterStride<>> SkipMap;
  SkipMap tmp(buffer, 0, 0, Eigen::OuterStride<>(1));
  for(auto Xsh : shell_range(Xblk)) {
    new (&tmp) SkipMap(buffer + block_offset,
        ish.nbf*jsh.nbf, Xsh.nbf, Eigen::OuterStride<>(Xblk.nbf)
    );
    ints_to_eigen_map(
        ish, jsh, Xsh,
        ints, int_type,
        tmp
    );
    block_offset += Xsh.nbf;
  }
  return rv;
}

template <typename MapType>
void
CADFCLHF::ints_to_eigen_map(
    const ShellData& ish,
    const ShellData& jsh,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type,
    MapType& out_map
){
  int block_offset = 0;
  const Eigen::Map<const RowMatrix> buffmap(ints->buffer(int_type), ish.nbf*jsh.nbf, Xsh.nbf);
  ints->compute_shell(ish, jsh, Xsh);
  out_map = buffmap;
  ints_computed_locally_ += ish.nbf * jsh.nbf * Xsh.nbf;
}


template <typename ShellRange>
Eigen::Map<CADFCLHF::ThreeCenterIntContainer>
CADFCLHF::ints_to_eigen_map(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type,
    double* buffer
){
  Eigen::Map<ThreeCenterIntContainer> rv(
      buffer,
      iblk.nbf * jsh.nbf, Xsh.nbf
  );
  int block_offset = 0;
  for(auto ish : shell_range(iblk)) {
    ints_to_buffer(
        ish, jsh, Xsh,
        ish.nbf, jsh.nbf, Xsh.nbf,
        ints, int_type,
        buffer + block_offset * jsh.nbf * Xsh.nbf
    );
    block_offset += ish.nbf;
  }
  return rv;
}

template <typename ShellRange1, typename ShellRange2>
Eigen::Map<CADFCLHF::ThreeCenterIntContainer>
CADFCLHF::ints_to_eigen_map(
    const ShellBlockData<ShellRange1>& iblk,
    const ShellData& jsh,
    const ShellBlockData<ShellRange2>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type,
    double* buffer
){
  Eigen::Map<ThreeCenterIntContainer> rv(
      buffer,
      iblk.nbf * jsh.nbf, Xblk.nbf
  );
  if(Xblk.nshell == 1) {
    const auto& Xsh = Xblk.first_shell;
    int block_offset = 0;
    for(auto&& ish : shell_range(iblk)) {
      ints_to_buffer(
          ish, jsh, Xsh,
          ish.nbf, jsh.nbf, Xsh.nbf,
          ints, int_type,
          buffer + block_offset * jsh.nbf * Xsh.nbf
      );
      block_offset += ish.nbf;
    }
    return rv;
  }
  else {
    int block_offset = 0;
    for(auto&& ish : shell_range(iblk)) {
      int Xblk_offset = 0;
      for(auto&& Xsh : shell_range(Xblk)) {
        ints_to_buffer(
            ish, jsh, Xsh,
            ish.nbf, jsh.nbf, Xsh.nbf,
            ints, int_type,
            buffer + block_offset * jsh.nbf * Xblk.nbf + Xblk_offset,
            // stride (from first element of one row to first element of next)
            Xblk.nbf
        );
        Xblk_offset += Xsh.nbf;
      }
      block_offset += ish.nbf;
    }
    return rv;
  }
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellData& ish,
    const ShellBlockData<ShellRange>& jblk,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      ish.nbf * jblk.nbf, Xsh.nbf
  );
  int block_offset = 0;
  for(auto jsh : shell_range(jblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
    for(auto mu : function_range(ish)){
      rv->block(
          mu.bfoff_in_shell*jblk.nbf + block_offset, 0,
          jsh.nbf, Xsh.nbf
      ) = ints_ptr->block(mu.bfoff_in_shell*jsh.nbf, 0, jsh.nbf, Xsh.nbf);
    }
    block_offset += jsh.nbf;
  }
  return rv;
}




} // end namespace sc

#endif /* _chemistry_qc_scf_get_ints_impl_h */
