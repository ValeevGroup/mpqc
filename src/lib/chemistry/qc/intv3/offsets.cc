//
// offsets.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#include <stdlib.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>

using namespace sc;

/* Compute the shell offset. */
static int
shell_offset(Ref<GaussianBasisSet> cs, int off)
{
  // unit shells have null cs's
  if (cs.null()) return off + 1;
  return off + cs->nshell();
}

/* Compute the prim offset. */
static int
prim_offset(Ref<GaussianBasisSet> cs, int off)
{
  // unit shells have null cs's
  if (cs.null()) return off + 1;
  return off + cs->nprimitive();
}

/* Compute the func offset. */
static int
func_offset(Ref<GaussianBasisSet> cs, int off)
{
  // unit shells have null cs's
  if (cs.null()) return off + 1;
  return off + cs->nbasis();
}

/////////////////////////////////////////////////////////////////////////

/* This initializes the offset arrays for one electron integrals. */
void
Int1eV3::int_initialize_offsets1()
{
  int shell_offset1;
  int prim_offset1;
  int func_offset1;

  /* Shell offset arrays. */
  bs1_shell_offset_ = 0;
  shell_offset1 = shell_offset(bs1_,0);
  if (bs2_ != bs1_) {
      shell_offset(bs2_,shell_offset1);
      bs2_shell_offset_ = shell_offset1;
    }
  else {
      bs2_shell_offset_ = bs1_shell_offset_;
    }

  /* Prim offset arrays. */
  bs1_prim_offset_ = 0;
  prim_offset1 = prim_offset(bs1_,0);
  if (bs2_ != bs1_) {
      prim_offset(bs2_,prim_offset1);
      bs2_prim_offset_ = prim_offset1;
    }
  else {
      bs2_prim_offset_ = bs1_prim_offset_;
    }

  /* Func offset arrays. */
  bs1_func_offset_ = 0;
  func_offset1 = func_offset(bs1_,0);
  if (bs2_ != bs1_) {
      func_offset(bs2_,func_offset1);
      bs2_func_offset_ = func_offset1;
    }
  else {
      bs2_func_offset_ = bs1_func_offset_;
    }
  }

/* This is called to free the offsets. */
void
Int1eV3::int_done_offsets1()
{
}

/* Initialize the offset arrays for two electron integrals. */
void
Int2eV3::int_initialize_offsets2()
{
  int shell_offset1;
  int shell_offset2;
  int shell_offset3;
  int prim_offset1;
  int prim_offset2;
  int prim_offset3;
  int func_offset1;
  int func_offset2;
  int func_offset3;

  /* Shell offset arrays. */
  bs1_shell_offset_ = 0;

  shell_offset1 = shell_offset(bs1_,0);
  if (bs2_ == bs1_) {
      shell_offset2 = shell_offset1;
      bs2_shell_offset_ = bs1_shell_offset_;
    }
  else {
      shell_offset2 = shell_offset(bs2_,shell_offset1);
      bs2_shell_offset_ = shell_offset1;
    }

  if (bs3_ == bs1_) {
      shell_offset3 = shell_offset2;
      bs3_shell_offset_ = bs1_shell_offset_;
    }
  else if (bs3_ == bs2_) {
      shell_offset3 = shell_offset2;
      bs3_shell_offset_ = bs2_shell_offset_;
    }
  else {
      shell_offset3 = shell_offset(bs3_,shell_offset2);
      bs3_shell_offset_ = shell_offset2;
    }

  if (bs4_ == bs1_) {
      bs4_shell_offset_ = bs1_shell_offset_;
    }
  else if (bs4_ == bs2_) {
      bs4_shell_offset_ = bs2_shell_offset_;
    }
  else if (bs4_ == bs3_) {
      bs4_shell_offset_ = bs3_shell_offset_;
    }
  else {
      bs4_shell_offset_ = shell_offset3;
    }

  /* Prim offset arrays. */
  bs1_prim_offset_ = 0;

  prim_offset1 = prim_offset(bs1_,0);
  if (bs2_ == bs1_) {
      prim_offset2 = prim_offset1;
      bs2_prim_offset_ = bs1_prim_offset_;
    }
  else {
      prim_offset2 = prim_offset(bs2_,prim_offset1);
      bs2_prim_offset_ = prim_offset1;
    }

  if (bs3_ == bs1_) {
      prim_offset3 = prim_offset2;
      bs3_prim_offset_ = bs1_prim_offset_;
    }
  else if (bs3_ == bs2_) {
      prim_offset3 = prim_offset2;
      bs3_prim_offset_ = bs2_prim_offset_;
    }
  else {
      prim_offset3 = prim_offset(bs3_,prim_offset2);
      bs3_prim_offset_ = prim_offset2;
    }

  if (bs4_ == bs1_) {
      bs4_prim_offset_ = bs1_prim_offset_;
    }
  else if (bs4_ == bs2_) {
      bs4_prim_offset_ = bs2_prim_offset_;
    }
  else if (bs4_ == bs3_) {
      bs4_prim_offset_ = bs3_prim_offset_;
    }
  else {
      bs4_prim_offset_ = prim_offset3;
    }

  /* Func offset arrays. */
  bs1_func_offset_ = 0;

  func_offset1 = func_offset(bs1_,0);
  if (bs2_ == bs1_) {
      func_offset2 = func_offset1;
      bs2_func_offset_ = bs1_func_offset_;
    }
  else {
      func_offset2 = func_offset(bs2_,func_offset1);
      bs2_func_offset_ = func_offset1;
    }

  if (bs3_ == bs1_) {
      func_offset3 = func_offset2;
      bs3_func_offset_ = bs1_func_offset_;
    }
  else if (bs3_ == bs2_) {
      func_offset3 = func_offset2;
      bs3_func_offset_ = bs2_func_offset_;
    }
  else {
      func_offset3 = func_offset(bs3_,func_offset2);
      bs3_func_offset_ = func_offset2;
    }

  if (bs4_ == bs1_) {
      bs4_func_offset_ = bs1_func_offset_;
    }
  else if (bs4_ == bs2_) {
      bs4_func_offset_ = bs2_func_offset_;
    }
  else if (bs4_ == bs3_) {
      bs4_func_offset_ = bs3_func_offset_;
    }
  else {
      bs4_func_offset_ = func_offset3;
    }
  }

/* This is called to free the offsets. */
void
Int2eV3::int_done_offsets2()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
