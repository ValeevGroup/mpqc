//
// state_debug.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _libqc_state_debug_h
#define _libqc_state_debug_h

#define SavableState_DECLARE_debug(classname) \
  public: \
    static classname* restore_state(StateIn&); \
    void save_parent_state(StateOut&); \
    void restore_parent_state(StateIn&); \
  private:

#define SavableState_IMPL_debug(classname) \
    classname* classname::restore_state(StateIn&si) \
      { \
      ClassDesc* cd; \
      int version; \
      DescribedClass* dc; \
      classname* ss = 0; \
      int objnum = si.getpointer((void**)&dc); \
      cout << "objnum = " << objnum << ", pointer = 0x" \
           << setbase(16) << (void*)dc << endl;\
      if (objnum) { \
        si.get(&cd,version); \
        if (!cd) { \
          cerr << #classname << "::restore_state: cannot restore" << endl; \
          exit(1); \
          } \
        dc = cd->create_described_class(); \
        if (!dc) { \
          cerr << #classname << "::restore_state: create failed" << endl; \
          exit(1); \
          } \
        ss = classname::castdown(dc); \
        if (!ss) { \
          cerr << #classname << "::restore_state: castdown failed(1)" << endl;\
          exit(1); \
          } \
        si.havepointer(objnum,(void*)dc); \
        ss->restore_vbase_state(version,si); \
        ss->restore_data_state(version,si); \
        } \
      else if (dc) { \
        ss = classname::castdown(dc); \
        if (!ss) { \
          cerr << #classname << "::restore_state: castdown failed(2)" << endl;\
          exit(1); \
          } \
        } \
      return ss; \
      } \
    void classname::save_parent_state(StateOut&so) \
      { \
      so.put(classname::class_desc()->version()); \
      classname::save_data_state(so); \
      } \
    void classname::restore_parent_state(StateIn&si) \
      { \
      int version; \
      si.get(version); \
      classname::restore_data_state(version,si); \
      }

#endif
