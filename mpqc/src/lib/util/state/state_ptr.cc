//
// state_ptr.cc
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

#include <ctype.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/stateptrImplSet.h>
#include <util/state/statenumImplSet.h>
#include <util/state/classdImplMap.h>

int
StateIn::dir_getobject(RefSavableState &p, const char *name)
{
  int r=0;

  p = 0;

  if (!use_directory()) {
      cerr << "ERROR: StateIn: no directory to get object from" << endl;
      abort();
    }

  // find the class name and/or object number
  const char *colon = strrchr(name,':');
  int number = 1;
  char *classname = 0;
  if (colon == 0) {
      if (isdigit(*name)) number = atoi(name);
      else classname = strcpy(new char[strlen(name)+1], name);
    }
  else {
      number = atoi(&colon[1]);
      classname = strcpy(new char[strlen(name)+1], name);
      *strrchr(classname,':') = '\0';
    }

  const ClassDesc *cd = 0;
  if (classname) {
      cd = ClassDesc::name_to_class_desc(classname);
      if (!cd) {
          cerr << "ERROR: StateIn: class " << classname << " unknown" << endl;
          abort();
        }
      delete[] classname;
    }

  int classid;
  if (cd) classid = classidmap_->operator[]((ClassDesc*)cd);
  else classid = -1;

  Pix i;
  int nfound = 0;
  for (i=ps_->first(); i; ps_->next(i)) {
      if (classid == -1 || ps_->operator()(i).type == classid) nfound++;
      if (nfound == number) {
          if (ps_->operator()(i).ptr.nonnull()) {
              p = ps_->operator()(i).ptr;
            }
          else {
              seek(ps_->operator()(i).offset);
              r += getobject(p);
            }
          return r;
        }
    }

  return r;
}

int
StateIn::getobject(RefSavableState &p)
{
  int r=0;
  int refnum;
  r += get(refnum);
  if (refnum == 0) {
      // reference to null
      p = 0;
    }
  else {
      StateDataNum num(refnum);
      Pix ind = ps_->seek(num);
      if (ind == 0 && use_directory()) {
          cerr << "ERROR: StateIn: directory missing object number "
               << refnum << endl;
          abort();
        }
      if (ind == 0 || ps_->operator()(ind).ptr.null()) {
          // object has not yet been read in
          if (use_directory()) {
              seek(ps_->operator()(ind).offset);
              int trefnum;
              get(trefnum);
              if (trefnum != refnum) {
                  cerr << "StateIn: didn't find expected reference" << endl;
                  abort();
                }
            }
          const ClassDesc *cd;
          r += get(&cd);
          have_classdesc();
          nextobject(refnum);
          DescribedClass *dc = cd->create(*this);
          p = SavableState::castdown(dc);
          if (use_directory()) {
              ps_->operator()(ind).ptr = p;
            }
        }
      else {
          // object already exists
          p = ps_->operator()(ind).ptr;
          if (use_directory() && tell() == ps_->operator()(ind).offset) {
              seek(tell() + ps_->operator()(ind).size);
            }
        }
    }
  return r;
}

void
StateIn::nextobject(int objnum)
{
  expected_object_num_ = objnum;
}

void
StateIn::haveobject(const RefSavableState &p)
{
  if (expected_object_num_) {
      haveobject(expected_object_num_,p);
      expected_object_num_ = 0;
    }
}

void
StateIn::haveobject(int objnum,const RefSavableState &p)
{
  StateDataNum num(objnum,p);
  Pix ind = ps_->seek(num);
  if (ind == 0) {
      ps_->add(num);
    }
  else {
      ps_->operator()(ind).ptr = p;
    }
}

int
StateOut::putobject(const RefSavableState &p)
{
  int r=0;
  if (p.null()) {
      // reference to null
      r += put(0);
    }
  else {
      StateDataPtr dp(p);
      Pix ind = ps_->seek(dp);
      if (ind == 0 || copy_references_) {
          // object has not been written yet
          dp.num = next_object_number_++;
          dp.offset = tell();
          r += put(dp.num);
          const ClassDesc *cd = p->class_desc();
          r += put(cd);
          dp.type = classidmap_->operator[]((ClassDesc*)cd);
          if (!copy_references_) ps_->add(dp);
          have_classdesc();
          p->save_vbase_state(*this);
          p->save_data_state(*this);
          if (!copy_references_) {
              Pix ind = ps_->seek(dp);
              ps_->operator()(ind).size = tell() - ps_->operator()(ind).offset;
            }
        }
      else {
          // object has already been written
          r += put(ps_->operator()(ind).num);
        }
    }
  return r;
}

/////////////////////////////////////////////////////////////////////////////
// StateData members

StateData::StateData(int n): num(n), ptr(0)
{
  init();
}

StateData::StateData(const RefSavableState &p): num(0), ptr(p)
{
  init();
}

StateData::StateData(int n, const RefSavableState &p): num(n), ptr(p)
{
  init();
}

void
StateData::init()
{
  size = 0;
  type = 0;
  offset = 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
