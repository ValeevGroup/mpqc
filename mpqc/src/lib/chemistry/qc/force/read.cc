
extern "C" {
#include <stdio.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
}
#include <util/keyval/keyval.h>
#include <util/keyval/ipv2c.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "mpcsscf.gbl"
#include "mposscf.gbl"
#include "read.h"

void
dmt_force_osscf_keyval_init(KeyVal*keyval,FILE*fp)
{
  dmt_force_osscf_init1(keyval,fp);
}

void
dmt_force_csscf_keyval_init(KeyVal*keyval,FILE*fp)
{
  dmt_force_csscf_init1(keyval,fp);
}

int
dmt_force_read_and_bcast_boolean(KeyVal*keyval,FILE*fp,char*name,int*boolval)
{
  int errcod=0;
  if (mynode0()==0) {
      if (keyval) {
          int saved_boolval = *boolval;
          *boolval = keyval->booleanvalue(name);
          errcod = (keyval->error() != KeyVal::OK);
          if (errcod) {
              *boolval = saved_boolval;
            }
        }
      else {
          errcod = ip_boolean(name,boolval,0);
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
      if (keyval) {
          int saved_intval = *intval;
          *intval = keyval->intvalue(name);
          errcod = (keyval->error() != KeyVal::OK);
          if (errcod) {
              *intval = saved_intval;
            }
        }
      else {
          errcod = ip_data(name,"%d",intval,0);
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
      if (keyval) {
          char* saved_stringptr = *stringptr;
          *stringptr = keyval.pcharvalue(name);
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
      else {
          errcod = ip_string(name,stringptr,0);
        }
    }
  return errcod;
}

