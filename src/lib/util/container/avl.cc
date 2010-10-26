//
// avl.cc
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

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif
#include <util/container/avlset.h>

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
#define INST_COMP(T) \
  template int compare(const T &, const T &)

INST_COMP(int);
INST_COMP(long);
INST_COMP(double);
INST_COMP(char);
INST_COMP(unsigned char);

template class EAVLMMapNode<int, AVLMapNode<int, int> >;
template class EAVLMMap<int, AVLMapNode<int, int> >;
template class AVLMapNode<int, int>;
template class AVLMap<int, int>;
template class AVLSet<int>;

#endif
