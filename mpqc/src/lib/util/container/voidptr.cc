//
// voidptr.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/container/voidptr.h>

VoidPtr::VoidPtr(): ptr_(0) {}
VoidPtr::VoidPtr(void* d): ptr_(d) {}
VoidPtr::VoidPtr(const VoidPtr& d): ptr_(d.ptr_) {}
VoidPtr::~VoidPtr() {}

void* VoidPtr::getptr() { return ptr_; }
int VoidPtr::operator==(VoidPtr d) { return ptr_ == d.ptr_; }
int VoidPtr::operator<=(VoidPtr d) { return ptr_ <= d.ptr_; }
VoidPtr& VoidPtr::operator=(void* d) { ptr_ = d; return *this; }
VoidPtr& VoidPtr::operator=(VoidPtr d) { ptr_ = d.ptr_; return *this; }
VoidPtr::operator void*() { return ptr_; }
