//
// eri.cc
//
// Copyright (C) 2001 Edward Valeev
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

#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/eri.h>
#ifdef DMALLOC
#include <dmalloc.h>
#endif

#if LIBINT2_SUPPORT_ERI

#define STORE_PAIR_DATA 1

using namespace std;
using namespace sc;

EriLibint2::EriLibint2(Integral *integral,
		   const Ref<GaussianBasisSet>& b1,
		   const Ref<GaussianBasisSet>& b2,
		   const Ref<GaussianBasisSet>& b3,
		   const Ref<GaussianBasisSet>& b4,
		   size_t storage) :
  Int2eLibint2(integral,b1,b2,b3,b4,storage),
  Fm_Eval_(b1->max_angular_momentum() +
          b2->max_angular_momentum() +
          b3->max_angular_momentum() +
          b4->max_angular_momentum())
{
  // The static part of Libint's interface is automatically initialized in libint.cc
  int l1 = bs1_->max_angular_momentum();
  int l2 = bs2_->max_angular_momentum();
  int l3 = bs3_->max_angular_momentum();
  int l4 = bs4_->max_angular_momentum();
  int lmax = max(max(l1,l2),max(l3,l4));
  if (lmax > LIBINT2_MAX_AM_ERI) {
    throw LimitExceeded<int>("EriLibint2::EriLibint2() -- maxam of the basis is too high,\
 not supported by this libint2 library. Recompile libint2.",__FILE__,__LINE__,LIBINT2_MAX_AM_ERI,lmax);
  }

  /*--- Initialize storage ---*/
  const int max_num_prim_comb = bs1_->max_nprimitive_in_shell()*
    bs2_->max_nprimitive_in_shell()*
    bs3_->max_nprimitive_in_shell()*
    bs4_->max_nprimitive_in_shell();
  // need one Libint_t object for each primitive combination
  // if Libint2 does not support contractions, just allocate 1
#if LIBINT2_CONTRACTED_INTS
  Libint_.resize(max_num_prim_comb);
#else
  Libint_.resize(1);
#endif
  ConsumableResources::get_default_instance()->consume_memory(Libint_.size() * sizeof(Libint_[0]));

  const int max_cart_target_size = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell()*
    bs3_->max_ncartesian_in_shell()*bs4_->max_ncartesian_in_shell();
  const int max_target_size = bs1_->max_nfunction_in_shell()*bs2_->max_nfunction_in_shell()*
    bs3_->max_nfunction_in_shell()*bs4_->max_nfunction_in_shell();

  size_t storage_needed = LIBINT2_PREFIXED_NAME(libint2_need_memory_eri)(lmax) * sizeof(LIBINT2_REALTYPE);
  LIBINT2_PREFIXED_NAME(libint2_init_eri)(&Libint_[0],lmax,0);  // only need to initialize stack of the first Libint_t object
  manage_array(Libint_[0].stack, storage_needed/sizeof(LIBINT2_REALTYPE));

  target_ints_buffer_ = allocate<double>(max_target_size);
  cart_ints_ = allocate<double>(max_cart_target_size);
  if (bs1_->has_pure() || bs2_->has_pure() || bs3_->has_pure() || bs4_->has_pure() ||
      bs1_->max_ncontraction() != 1 || bs2_->max_ncontraction() != 1 ||
      bs3_->max_ncontraction() != 1 || bs4_->max_ncontraction() != 1) {
    sphharm_ints_ = allocate<double>(max_target_size);
    storage_needed += max_target_size*sizeof(double);
  }
  else {
    sphharm_ints_ = 0;
  }
  if (l1 || l2 || l3 || l4) {
    perm_ints_ = allocate<double>(max_target_size);
    storage_needed += max_target_size*sizeof(double);
  }
  else
    perm_ints_ = 0;

  // See if can store primitive-pair data
  size_t primitive_pair_storage_estimate = (bs1_->nprimitive()*bs2_->nprimitive() + 
    bs3_->nprimitive()*bs4_->nprimitive())*sizeof(prim_pair_t);
  //  ExEnv::errn() << scprintf("need %d bytes to store primitive pair data\n",primitive_pair_storage_estimate);
#if STORE_PAIR_DATA
  shell_pairs12_ = new ShellPairsLibint2(bs1_,bs2_);
  if ( (bs1_ == bs3_ && bs2_ == bs4_) /*||
       // if this is (ab|ba) case -- should i try to save storage?
       (bs1_ == bs4_ && bs2_ == bs3_)*/ )
    shell_pairs34_ = new ShellPairsLibint2(shell_pairs12_);
  else
    shell_pairs34_ = new ShellPairsLibint2(bs3_,bs4_);
  storage_needed += primitive_pair_storage_estimate;
#endif

  storage_used_ = storage_needed;
  // Check if storage_ > storage_needed
  check_storage_();

  int mmax = bs1_->max_angular_momentum() +
    bs2_->max_angular_momentum() +
    bs3_->max_angular_momentum() +
    bs4_->max_angular_momentum();
  Fm_table_ = new double[mmax+1];
}


EriLibint2::~EriLibint2()
{ 
  unmanage_array(Libint_[0].stack);
  LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&Libint_[0]);
  Libint_[0].stack = 0;
  ConsumableResources::get_default_instance()->release_memory(Libint_.size() * sizeof(Libint_[0]));
  deallocate(target_ints_buffer_);
  deallocate(cart_ints_);
  if (sphharm_ints_)
    deallocate(sphharm_ints_);
  if (perm_ints_)
    deallocate(perm_ints_);
#ifdef DMALLOC
  dmalloc_shutdown();
#endif
  delete[] Fm_table_;
}

size_t
EriLibint2::storage_required(const Ref<GaussianBasisSet>& b1,
			   const Ref<GaussianBasisSet>& b2,
			   const Ref<GaussianBasisSet>& b3,
			   const Ref<GaussianBasisSet>& b4)
{
  Ref<GaussianBasisSet> bs1 = b1;
  Ref<GaussianBasisSet> bs2 = b2;
  Ref<GaussianBasisSet> bs3 = b3;
  Ref<GaussianBasisSet> bs4 = b4;

  if (bs2 == 0)
    bs2 = bs1;
  if (bs3 == 0)
    bs3 = bs1;
  if (bs4 == 0)
    bs4 = bs1;

  int l1 = bs1->max_angular_momentum();
  int l2 = bs2->max_angular_momentum();
  int l3 = bs3->max_angular_momentum();
  int l4 = bs4->max_angular_momentum();
  int lmax = max(max(l1,l2),max(l3,l4));

  size_t storage_required = storage_required_(bs1,bs2,bs3,bs4);

  const int max_num_prim_comb = bs1->max_nprimitive_in_shell()*
    bs2->max_nprimitive_in_shell()*
    bs3->max_nprimitive_in_shell()*
    bs4->max_nprimitive_in_shell();
#if LIBINT2_CONTRACTED_INTS
  storage_required += max_num_prim_comb * sizeof(Libint_t);
#else
  storage_required += sizeof(Libint_t);
#endif

  const int max_cart_target_size = bs1->max_ncartesian_in_shell()*bs2->max_ncartesian_in_shell()*
    bs3->max_ncartesian_in_shell()*bs4->max_ncartesian_in_shell();
  const int max_target_size = bs1->max_nfunction_in_shell()*bs2->max_nfunction_in_shell()*
    bs3->max_nfunction_in_shell()*bs4->max_nfunction_in_shell();

  storage_required += LIBINT2_PREFIXED_NAME(libint2_need_memory_eri)(lmax) * sizeof(LIBINT2_REALTYPE);

  if (bs1->has_pure() || bs2->has_pure() || bs3->has_pure() || bs4->has_pure() ||
      bs1->max_ncontraction() != 1 || bs2->max_ncontraction() != 1 ||
      bs3->max_ncontraction() != 1 || bs4->max_ncontraction() != 1) {
    storage_required += max_target_size*sizeof(double);
  }

  if (l1 || l2 || l3 || l4) {
    storage_required += max_target_size*sizeof(double);
  }

  // See if can store primitive-pair data
  size_t primitive_pair_storage_estimate = (bs1->nprimitive()*bs2->nprimitive() + 
    bs3->nprimitive()*bs4->nprimitive())*sizeof(prim_pair_t);
#if STORE_PAIR_DATA
  storage_required += primitive_pair_storage_estimate;
#endif

  return storage_required;
}

#endif // if LIBINT2_SUPPORT_ERI

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
