//
// messpvm.cc
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

#include <pvm3.h>
// extern "C" {
//     int pvm_catchout(FILE*fp);
// }
#include <util/keyval/keyval.h>
#include <util/group/messpvm.h>

#define CLASSNAME PVMMessageGrp
#define PARENTS public MessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
PVMMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

PVMMessageGrp::PVMMessageGrp(const RefKeyVal&keyval):
  MessageGrp(keyval)
{
  int i;

  int master_is_node = keyval->booleanvalue("master_is_node");
  if (keyval->error() != KeyVal::OK) master_is_node = 1;
  if (master_is_node != 0) master_is_node = 1;

  int nprocs = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) nprocs = 1;

  char* task = keyval->pcharvalue("executable");
  if (keyval->error() != KeyVal::OK) {
      ExEnv::err() << scprintf("PVMMessageGrp(KeyVal): \"executable\" not given\n");
      abort();
    }

  int nhosts = keyval->count("hosts");
  char **hosts;
  if (keyval->error() == KeyVal::OK && nhosts > 0) {
      hosts = new char*[nhosts];
      for (i=0; i<nhosts; i++) {
          hosts[i] = keyval->pcharvalue("hosts",i);
        }
    }
  else {
      hosts = 0;
    }

  char *envstr = "MessageGrp=PVMMessageGrp";
  char *tmpstr = (char*)malloc(strlen(envstr)+1);
  strcpy(tmpstr,envstr);
  putenv(tmpstr);

  envstr = "PVM_EXPORT=MessageGrp";
  tmpstr = (char*)malloc(strlen(envstr)+1);
  strcpy(tmpstr,envstr);
  putenv(tmpstr);

  int numt = 0;
  int n_to_spawn = (master_is_node?nprocs - 1: nprocs);
  char **argv = 0;
  tids = new int[nprocs];
  // not sure about this under c++
  //pvm_catchout(stdout);
  if (master_is_node) {
      tids[0] = pvm_mytid();
    }
  if (hosts) {
      int ihost = master_is_node;
      for (i=master_is_node; i<nprocs; i++) {
          if (ihost >= nhosts) ihost = 0;
          numt += pvm_spawn(task, argv,
                            PvmTaskHost, hosts[ihost], 1, &tids[i]);
          ihost++;
        }
    }
  else {
      numt += pvm_spawn(task, argv,
                        PvmTaskDefault, 0,
                        n_to_spawn, &tids[master_is_node]);
    }

  if (numt < n_to_spawn) {
      ExEnv::err() << scprintf("PVMMessageGrp(KeyVal): failed to spawn all processes\n");
      ExEnv::err() << scprintf(" numt = %d, n_to_spawn = %d\n", numt, n_to_spawn);
      for (i=0; i<nprocs; i++) {
          ExEnv::err() << scprintf("tids[%d] = %d\n", i, tids[i]);
        }
      abort();
    }

  for (i=0; i<nhosts; i++) {
      delete[] hosts[i];
    }
  if (hosts) delete[] hosts;

  // send the tids to all of the spawned processes
  for (i=master_is_node; i<nprocs; i++) {
      pvm_psend(tids[i], 2, (char*) &nprocs, 1, PVM_INT);
      pvm_psend(tids[i], 2, (char*) tids, nprocs, PVM_INT);
    }
  if (!master_is_node) {
      pvm_exit();
      exit(0);
    }

  initialize(0,nprocs);
}

PVMMessageGrp::PVMMessageGrp()
{
  int rtid, rtag, rlen;
  int nproc;
  int parent = pvm_parent();
  pvm_precv(parent, 2, (void*) &nproc, 1, PVM_INT, &rtid, &rtag, &rlen);
  tids = new int[nproc];
  pvm_precv(parent, 2, (void*) tids, nproc, PVM_INT, &rtid, &rtag, &rlen);

  int mytid = pvm_mytid();
  for (int me=0; me<nproc; me++) {
      if (tids[me] == mytid) break;
    }

  initialize(me, nproc);
}

PVMMessageGrp::~PVMMessageGrp()
{
  delete[] tids;
  pvm_exit();
}

void
PVMMessageGrp::raw_send(int target, void* data, int nbyte)
{
  pvm_psend(tids[target], 0, (char*)data, nbyte, PVM_BYTE);
}

void
PVMMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  pvm_precv(tids[sender], 0, (char*)data, nbyte, PVM_BYTE,
            &rtid, &rtag, &rlen);
}

void
PVMMessageGrp::raw_sendt(int target, int type, void* data, int nbyte)
{
  pvm_psend(tids[target], (type<<1) + 1, (char*)data, nbyte, PVM_BYTE);
}

void
PVMMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  pvm_precv(-1, (type<<1) + 1, (char*)data, nbyte, PVM_BYTE,
            &rtid, &rtag, &rlen);
}

int
PVMMessageGrp::probet(int type)
{
  int bufid = pvm_probe(-1, (type<<1) + 1);
  if (bufid > 0) {
      pvm_bufinfo(bufid, &rlen, &rtag, &rtid);
      return 1;
    }
  return 0;
}

int
PVMMessageGrp::last_source()
{
  // convert the tid to a node number
  for (int i=0; i<n(); i++) {
      if (tids[i] == rtid) break;
    }
  return i;
}

int
PVMMessageGrp::last_size()
{
  return rlen;
}

int
PVMMessageGrp::last_type()
{
  return rtag>>1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
