//
// tmplinst.cc
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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
/////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif
#include <util/container/avlmap.h>
#include <util/state/state.h>
#include <util/state/stateio.h>

using namespace sc;

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

template class EAVLMMapNode<int, AVLMapNode<int, StateInData> >;
template class EAVLMMap<int, AVLMapNode<int, StateInData> >;
template class AVLMapNode<int, StateInData>;
template class AVLMap<int, StateInData>;

template class EAVLMMapNode<Ref<SavableState>, AVLMapNode<Ref<SavableState>, StateOutData> >;
template class EAVLMMap<Ref<SavableState>, AVLMapNode<Ref<SavableState>, StateOutData> >;
template class AVLMapNode<Ref<SavableState>, StateOutData>;
template class AVLMap<Ref<SavableState>, StateOutData>;

template class EAVLMMapNode<ClassDescP, AVLMapNode<ClassDescP,int> >;
template class EAVLMMap<ClassDescP, AVLMapNode<ClassDescP,int> >;
template class AVLMapNode<ClassDescP,int>;
template class AVLMap<ClassDescP,int>;

template class EAVLMMapNode<int, AVLMapNode<int,StateClassData> >;
template class EAVLMMap<int, AVLMapNode<int,StateClassData> >;
template class AVLMapNode<int,StateClassData>;
template class AVLMap<int,StateClassData>;

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
