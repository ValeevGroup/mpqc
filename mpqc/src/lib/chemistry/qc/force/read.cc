//
// read.cc
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

extern "C" {
#include <stdio.h>
}
#include <util/group/picl.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/force/mpcsscf.gbl>
#include <chemistry/qc/force/mposscf.gbl>
#include <chemistry/qc/force/read.h>

void
dmt_force_osscf_keyval_init(KeyVal*keyval,FILE*fp)
{
  if (!keyval && mynode0() == 0) {
      fprintf(stderr, "dmt_force_osscf_keyval_init got null keyval\n");
      abort();
    }
  dmt_force_osscf_init1(keyval,fp);
}

void
dmt_force_csscf_keyval_init(KeyVal*keyval,FILE*fp)
{
  if (!keyval && mynode0() == 0) {
      fprintf(stderr, "dmt_force_csscf_keyval_init got null keyval\n");
      abort();
    }
  dmt_force_csscf_init1(keyval,fp);
}

int
dmt_force_read_and_bcast_boolean(KeyVal*keyval,FILE*fp,char*name,int*boolval)
{
  int errcod=0;
  if (mynode0()==0) {
      int saved_boolval = *boolval;
      *boolval = keyval->booleanvalue(name);
      errcod = (keyval->error() != KeyVal::OK);
      if (errcod) {
          *boolval = saved_boolval;
        }
    }
  bcast0(boolval,sizeof(int),mtype_get(),0);
  if (mynode0()==0) fprintf(fp,"  :force:%s = %d\n",name,*boolval);
  return errcod;
}

int
dmt_force_read_and_bcast_int(KeyVal*keyval,FILE*fp,char*name,int*intval)
{
  int errcod=0;
  if (mynode0()==0) {
      int saved_intval = *intval;
      *intval = keyval->intvalue(name);
      errcod = (keyval->error() != KeyVal::OK);
      if (errcod) {
          *intval = saved_intval;
        }
    }
  bcast0(intval,sizeof(int),mtype_get(),0);
  if (mynode0()==0) fprintf(fp,"  :force:%s = %d\n",name,*intval);
  return errcod;
}


int
dmt_force_read_string(KeyVal*keyval,FILE*fp,char*name,char**stringptr)
{
  int errcod=0;
  if (mynode0()==0) {
      char* saved_stringptr = *stringptr;
      *stringptr = keyval->pcharvalue(name);
      errcod = (keyval->error() != KeyVal::OK);
      if (errcod) {
          *stringptr = saved_stringptr;
        }
      else {
          // convert new'ed memory to malloc'ed memory
          char* tmp = *stringptr;
          *stringptr = (char*) malloc(strlen(tmp) + 1);
          strcpy(*stringptr,tmp);
          delete[] tmp;
        }
    }
  return errcod;
}

void
dmt_get_csscf_force(StateIn& s)
{
  int *opts=0;
  if (mynode0()==0)
    s.get(opts);
  else
    opts = new int[17];
  bcast0(opts,sizeof(int)*17,mtype_get(),0);
  
  dmt_force_csscf_put_options(opts);
  delete[] opts;
}

void
dmt_put_csscf_force(StateOut& s)
{
  if (mynode0()==0) {
    int *opts=new int[17];
    dmt_force_csscf_get_options(opts);
    s.put(opts,17);
    delete[] opts;
  }
}

void
dmt_get_osscf_force(StateIn& s)
{
  int *opts=0;
  if (mynode0()==0)
    s.get(opts);
  else
    opts = new int[17];
  bcast0(opts,sizeof(int)*17,mtype_get(),0);
  
  dmt_force_osscf_put_options(opts);
  delete[] opts;
}

void
dmt_put_osscf_force(StateOut& s)
{
  if (mynode0()==0) {
    int *opts=new int[17];
    dmt_force_osscf_get_options(opts);
    s.put(opts,17);
    delete[] opts;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
