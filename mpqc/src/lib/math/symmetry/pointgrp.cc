//
// pointgrp.cc
//
// Modifications are
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

/* pointgrp.cc -- implementation of the point group classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      June, 1993
 */

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <math/symmetry/pointgrp.h>

////////////////////////////////////////////////////////////////////////

#define CLASSNAME PointGroup
#define PARENTS public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PointGroup::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

PointGroup::PointGroup()
  : symb(0)
{
  set_symbol("c1");
  frame(0,0) = frame(1,1) = frame(2,2) = 1;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s)
  : symb(0)
{
  set_symbol(s);
  frame(0,0) = frame(1,1) = frame(2,2) = 1;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s, SymmetryOperation& so)
  : symb(0)
{
  set_symbol(s);
  frame = so;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s, SymmetryOperation& so, Point& or)
  : symb(0)
{
  set_symbol(s);
  frame = so;
  origin_ = or;
}

PointGroup::PointGroup(const RefKeyVal& kv)
  : symb(0)
{
  if (kv->exists("symmetry")) {
    char *tmp = kv->pcharvalue("symmetry");
    set_symbol(tmp);
    delete[] tmp;
  }
  else
    set_symbol("c1");

  if (kv->exists("symmetry_frame")) {
    for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++) 
        frame(i,j) = kv->doublevalue("symmetry_frame",i,j);
  } else {
    frame(0,0) = frame(1,1) = frame(2,2) = 1;
  }

  if (kv->exists("origin")) {
    for (int i=0; i < 3; i++)
      origin_[i] = kv->doublevalue("origin",i);
  } else {
    origin_[0] = origin_[1] = origin_[2] =0;
  }
}

PointGroup::PointGroup(StateIn& si) :
  symb(0),
  SavableState(si),
  origin_(si)
{
  si.getstring(symb);
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      si.get(frame(i,j));
}

PointGroup::PointGroup(const PointGroup& pg)
  : symb(0)
{
  *this = pg;
}

PointGroup::~PointGroup()
{
  if (symb) { delete[] symb; symb=0; }
}

PointGroup&
PointGroup::operator=(const PointGroup& pg)
{
  set_symbol(pg.symb);
  frame = pg.frame;
  origin_ = pg.origin_;
  return *this;
}

void
PointGroup::set_symbol(const char *sym)
{
  if (sym) {
    if (symb) delete[] symb;
    int len;
    symb = new char[(len=strlen(sym))+1];
    for (int i=0; i<len; i++) symb[i] = (char) tolower(sym[i]);
    symb[len] = '\0';
  } else {
    set_symbol("c1");
  }
}

void
PointGroup::save_data_state(StateOut& so)
{
  origin_.save_object_state(so);
  so.putstring(symb);

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      so.put(frame(i,j));
}

CharacterTable
PointGroup::char_table() const
{
  CharacterTable ret(symb,frame);
  return ret;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
