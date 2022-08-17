//
// messmpi.cc
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

#include <algorithm>

#include <stdio.h> // for sprintf
#include <unistd.h> // for fchdir etc.
#include <fcntl.h> // for open on AIX

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
extern int MPI_Initialized(int *); // used to be missing in old mpi.h

#include <util/keyval/keyval.h>
#include <util/group/messmpi.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>

MPI_Comm global_commgrp;

using namespace std;
using namespace sc;

// OP_COMMUTES is zero to work around a bug in MPI/Pro 1.5b5 and earlier
#define OP_COMMUTES 1

///////////////////////////////////////////////////////////////////////

int MPIMessageGrp::nmpi_grps=0;
ThreadLock* MPIMessageGrp::grplock=0;
bool MPIMessageGrp::mpi_init_called=false;

static
void
print_error_and_abort(int me, int mpierror)
{
  char msg[MPI_MAX_ERROR_STRING+1];
  int size;
  MPI_Error_string(mpierror, msg, &size);
  msg[size] = '\0';
  ExEnv::outn() << me << ": " << msg << endl;
  ExEnv::outn().flush();
  //MPI_Abort(MPI_COMM_WORLD, mpierror);
}

static
const char *
mpi_thread_string(int level)
{
  switch (level) {
#ifdef HAVE_MPI_INIT_THREAD
  case MPI_THREAD_SINGLE:
      return "MPI_THREAD_SINGLE";
  case MPI_THREAD_FUNNELED:
      return "MPI_THREAD_FUNNELED";
  case MPI_THREAD_SERIALIZED:
      return "MPI_THREAD_SERIALIZED";
  case MPI_THREAD_MULTIPLE:
      return "MPI_THREAD_MULTIPLE";
#endif
  default:
      return "unknown";
    }
}

///////////////////////////////////////////////////////////////////////
// The MPIMessageGrp class

static ClassDesc MPIMessageGrp_cd(
  typeid(MPIMessageGrp),"MPIMessageGrp",1,"public MessageGrp",
  create<MPIMessageGrp>, create<MPIMessageGrp>, 0);

MPIMessageGrp::MPIMessageGrp()
{
  init(MPI_COMM_WORLD);
}

MPIMessageGrp::MPIMessageGrp(MPI_Comm comm)
{
  init(comm);
}

MPIMessageGrp::MPIMessageGrp(int *argc, char ***argv)
{
  init(MPI_COMM_WORLD,argc,argv);
}

MPIMessageGrp::MPIMessageGrp(const Ref<KeyVal>& keyval):
  MessageGrp(keyval)
{
  if (keyval->exists("argv")) {
      int argc = keyval->count("argv");
      char **argv = new char*[argc+1];
      argv[argc] = 0;
      for (int arg=0; arg<argc; arg++) {
          std::string str = keyval->stringvalue("argv",arg);
          // this is a minor memory leak
          argv[arg] = strcpy(new char[str.size()+1], str.c_str());
        }
      init(MPI_COMM_WORLD, &argc, &argv);
    }
  else {
      init(MPI_COMM_WORLD);
    }

  if (keyval->booleanvalue("errors_return")) {
      if (me()==0)
          ExEnv::outn() << indent << "MPIMessageGrp: errors_return is true" << endl;
      MPI_Comm_set_errhandler(commgrp, MPI_ERRORS_RETURN);
    }

  if (debug_) {
      ExEnv::outn() << indent << "MPIMessageGrp: KeyVal CTOR: done" << endl;
    }
}

void
MPIMessageGrp::init(MPI_Comm comm, int *argc, char ***argv)
{
  int me, nproc;

  if (debug_) {
      ExEnv::outn() << "MPIMessageGrp::init: entered" << endl;
    }

  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
      int tmp_argc;
      char **tmp_argv;
      int *inits_argc;
      char ***inits_argv;
      if (argc && argv) {
          inits_argc = argc;
          inits_argv = argv;
        }
      else {
          tmp_argc = 0;
          tmp_argv = new char*[tmp_argc+1];
          tmp_argv[tmp_argc] = 0;
          inits_argc = &tmp_argc;
          inits_argv = &tmp_argv;
        }
      // This dot business is to work around problems with some MPI
      // implementations.
      int dot = open(".",O_RDONLY);
      if (debug_) {
          ExEnv::outn() << indent
               << "Calling MPI_Init with";
          for (int i=0; i<*inits_argc; i++) {
              ExEnv::outn() << " " << *inits_argv[i];
            }
          ExEnv::outn() << endl;
        }
#ifdef HAVE_MPI_INIT_THREAD
      int provided, desired = MPQC_MPI_THREAD_LEVEL;
      MPI_Init_thread(inits_argc, inits_argv, desired, &provided);
      int me;
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      if (provided != desired) {
          ExEnv::out0() << indent
                        << "WARNING: desired "
                        << mpi_thread_string(desired)
                        << " MPI threading support but got "
                        << mpi_thread_string(provided)
                        << endl;
        }
#else
      MPI_Init(inits_argc, inits_argv);
#endif
      mpi_init_called = true;
#ifdef HAVE_FCHDIR
      fchdir(dot);
#endif
      close(dot);
    }

  MPI_Comm_dup(comm, &commgrp);
  global_commgrp = commgrp;

  MPI_Comm_set_errhandler(commgrp, MPI_ERRORS_ARE_FATAL);

  if (!nmpi_grps) {
    threadgrp = ThreadGrp::get_default_threadgrp();
    Ref<ThreadLock> glock = threadgrp->new_lock();
    // to avoid auto-deletion of grplock
    glock->reference();
    grplock = glock.pointer();
  }
  {
    grplock->lock();
    nmpi_grps++;
    grplock->unlock();
  }

  MPI_Comm_rank(commgrp,&me);
  MPI_Comm_size(commgrp, &nproc);
  bufsize = 4000000;
  buf = 0;
  //buf = (void*) new char[bufsize];
  //MPI_Buffer_attach(buf,bufsize);

  if (getenv("MPIMESSAGEGRP_MESSAGEGRP_COLLECTIVES"))
      use_messagegrp_collectives_ = true;
  else
      use_messagegrp_collectives_ = false;
  
  initialize(me, nproc);

  //MPIL_Trace_on();

  if (debug_) {
      ExEnv::outn() << me << ": MPIMessageGrp::init: done" << endl;
    }

  SCFormIO::init_mp(me);
}

MPIMessageGrp::~MPIMessageGrp()
{
  //MPIL_Trace_off();
  //MPI_Buffer_detach(&buf, &bufsize);
  delete[] (char*) buf;

  {
    Ref<ThreadLock> lock = grplock;
    lock->lock();
    nmpi_grps--;
    if (nmpi_grps==0 && mpi_init_called) {
      MPI_Finalize();
      mpi_init_called = false;
      // grplock will be cleaned up after this
      grplock->dereference();
    }
    else
      MPI_Comm_free(&commgrp);
    lock->unlock();
  }
}

Ref<MessageGrp> MPIMessageGrp::clone(void)
{
  Ref<MessageGrp> mgrp = new MPIMessageGrp;
  return mgrp;
}

Ref<MessageGrp> MPIMessageGrp::split(int grpkey, int rankkey)
{
  Ref<MessageGrp> mgrp;
  MPI_Comm newcomm;
  int color = grpkey;
  if (color < 0) color = MPI_UNDEFINED;
  MPI_Comm_split(commgrp, color, rankkey, &newcomm);
  if (grpkey >=0 ) mgrp = new MPIMessageGrp(newcomm);
  return mgrp;
}

Ref<MessageGrp> MPIMessageGrp::subset(const std::set<int> &s)
{
  std::vector<int> v(s.size());
  std::copy(s.begin(), s.end(), v.begin());

  Ref<MessageGrp> mgrp;

  MPI_Group grp, grp2;
  MPI_Comm_group(commgrp, &grp);
  MPI_Group_incl(grp, v.size(), &v.front(), &grp2);

  MPI_Comm newcomm;
  MPI_Comm_create(commgrp, grp2, &newcomm);
  if (s.find(me()) != s.end()) mgrp = new MPIMessageGrp(newcomm);

  MPI_Group_free(&grp);
  MPI_Group_free(&grp2);

  return mgrp;
}

void
MPIMessageGrp::raw_send(int target, const void* data, int nbyte)
{
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Send"
                       "(0x%08x, %5d, MPI_BYTE, %3d, 0, commgrp)",
                       me(), data, nbyte, target)
           << endl;
    }
  int ret;
  ret = MPI_Send(const_cast<void*>(data),nbyte,MPI_BYTE,target,0,commgrp);
  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_send("
          << target << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) ExEnv::outn() << scprintf("%3d: sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recv(int sender, void* data, int nbyte,
                        MessageInfo *info)
{
  MPI_Status status;
  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Recv"
                       "(0x%08x, %5d, MPI_BYTE, %3d, 0, commgrp,)",
                       me(), data, nbyte, sender)
           << endl;
    }
  int ret;
  ret = MPI_Recv(data,nbyte,MPI_BYTE,sender,0,commgrp,&status);
  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_recv("
          << sender << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  set_sender(info,status.MPI_SOURCE);
  set_type(info,-1);
  set_nbyte(info,nbyte);
  if (debug_) {
      int rnode = status.MPI_SOURCE;
      int rlen = nbyte;
      ExEnv::outn() << scprintf("%3d: recvd %d bytes\n", me(), rlen) << endl;
    }
}

void
MPIMessageGrp::raw_sendt(int target, int type, const void* data, int nbyte,
                         bool rcvrdy)
{
  type = type*2 + 1;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Send"
                       "(0x%08x, %5d, MPI_BYTE, %3d, %5d, commgrp)",
                       me(), data, nbyte, target, type)
           << endl;
    }
  int ret;
  if (rcvrdy) {
      ret = MPI_Rsend(const_cast<void*>(data),nbyte,MPI_BYTE,target,type,commgrp);
    }
  else {
      ret = MPI_Send(const_cast<void*>(data),nbyte,MPI_BYTE,target,type,commgrp);
    }
  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_sendt("
          << target << "," << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) ExEnv::outn() << scprintf("%3d: sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recvt(int sender, int type, void* data, int nbyte,
                         MessageInfo *info)
{
  MPI_Status status;
  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Recv(0x%08x, %5d, MPI_BYTE, "
                       "%5d, %5d, commgrp,)",
                       me(), data, nbyte, sender, type)
           << endl;
    }
  int ret;
  ret = MPI_Recv(data,nbyte,MPI_BYTE,sender,type,commgrp,&status);
  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_recvt("
          << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (info) {
      set_sender(info,status.MPI_SOURCE);
      set_type(info,(status.MPI_TAG-1)/2);
      set_nbyte(info,nbyte);
    }
  if (debug_) {
      int rnode = status.MPI_SOURCE;
      int rtag = status.MPI_TAG;
      int rlen = nbyte;
      ExEnv::outn() << scprintf("%3d: recvd %d bytes from %d with tag %d\n",
                       me(), rlen, rnode, rtag) << endl;
    }
}

void
MPIMessageGrp::raw_nb_sendt(int target, int type, const void* data, int nbyte,
                            MessageHandle &mh,
                            bool rcvrdy)
{
  int tag = type*2 + 1;
  int ret;

  MessageHandleData *mhd = new MessageHandleData(nbyte);

  if (rcvrdy) {
      ret = MPI_Irsend(const_cast<void*>(data),
                       nbyte,MPI_BYTE,target,tag,commgrp,&mhd->req);
    }
  else {
      ret = MPI_Isend(const_cast<void*>(data),
                      nbyte,MPI_BYTE,target,tag,commgrp,&mhd->req);
    }

  if (debug_) {
      ExEnv::outn() << scprintf("%3d: isend %d bytes to %d with tag %d\n",
                       me(), nbyte, target, tag) << endl;
    }

  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_nb_sendt("
          << target << "," << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }

  set_id(&mh, mhd);
}

void
MPIMessageGrp::raw_nb_recvt(int sender, int type, void* data, int nbyte,
                            MessageHandle &mh)
{
  if (sender == -1) sender = MPI_ANY_SOURCE;

  int tag;
  if (type == -1) tag = MPI_ANY_TAG;
  else tag = (type<<1) + 1;

  MessageHandleData *mhd = new MessageHandleData(nbyte);

  int ret;
  ret = MPI_Irecv(data,nbyte,MPI_BYTE,sender,tag,commgrp,&mhd->req);

  if (debug_) {
      ExEnv::outn() << scprintf("%3d: irecv %d bytes from %d with tag %d\n",
                       me(), nbyte, sender, tag) << endl;
    }

  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_nb_recvt("
          << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }

  set_id(&mh, mhd);
}

void
MPIMessageGrp::wait(const MessageHandle&mh,
                    MessageInfo *info)
{
  MessageHandleData *mhd = static_cast<MessageHandleData*>(get_id(&mh));
  MPI_Status status;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: wait\n",
                       me()) << endl;
    }
  int ret = MPI_Wait(&mhd->req,&status);
  if (ret != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::wait: mpi error:"
                    << endl;
      print_error_and_abort(me(), ret);
    }
  if (info) {
      set_sender(info,status.MPI_SOURCE);
      set_type(info,(status.MPI_TAG-1)/2);
      set_nbyte(info,mhd->nbyte);
    }
}

int
MPIMessageGrp::probet(int sender, int type,
                      MessageInfo *info)
{
  int flag;
  MPI_Status status;

  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  int ret;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Iprobe(%5d, %5d, commgrp, "
                       "&flag, &status)", me(), sender, type)
           << endl;
    }
  if ((ret = MPI_Iprobe(sender,type,commgrp,&flag,&status))
      != MPI_SUCCESS ) {
      ExEnv::outn() << me() << ": MPIMessageGrp::probet("
          << type << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (flag) {
      int nbyte;
      MPI_Get_count(&status, MPI_BYTE, &nbyte);
      set_sender(info,status.MPI_SOURCE);
      set_type(info,(status.MPI_TAG-1)/2);
      set_nbyte(info,nbyte);
      return 1;
    }
  else {
      if (info) {
          set_sender(info,-1);
          set_type(info,-1);
          set_nbyte(info,0);
        }
    }
    
  return 0;
}

void
MPIMessageGrp::sync()
{
  int ret;
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Barrier(commgrp)", me()) << endl;
    }
  if ((ret = MPI_Barrier(commgrp)) != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::sync(): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
}

#define REDUCEMEMBER(name, type, mpitype) \
static GrpReduce<type>* name ## reduceobject; \
extern "C" void \
name ## reduce(void*b, void*a, int*len, MPI_Datatype*datatype) \
{ \
  name ## reduceobject->reduce((type*)a, (type*)b, *len); \
} \
void \
MPIMessageGrp::reduce(type*d, int n, GrpReduce<type>&r, \
                      type*scratch, int target) \
{ \
  if (use_messagegrp_collectives_) { \
      MessageGrp::reduce(d,n,r,scratch,target); \
      return; \
    } \
 \
  name ## reduceobject = &r; \
 \
  MPI_Op op; \
  MPI_Op_create(name ## reduce, OP_COMMUTES, &op); \
 \
  type *work; \
  if (!scratch) work = new type[n]; \
  else work = scratch; \
 \
  int ret; \
 \
  if (target == -1) { \
      if (debug_) { \
          ExEnv::outn() << scprintf("%3d: MPI_Allreduce" \
          "(0x%08x, 0x%08x, %5d, %3d, op, commgrp)", \
          me(), d, work, n, mpitype) \
               << endl; \
        } \
      ret = MPI_Allreduce(d, work, n, mpitype, op, commgrp); \
      if (debug_) \
        ExEnv::outn() << scprintf("%3d: done with Allreduce", me()) << endl; \
    } \
  else { \
      if (debug_) { \
          ExEnv::outn() << scprintf("%3d: MPI_Reduce" \
          "(0x%08x, 0x%08x, %5d, %3d, op, %3d, commgrp)", \
          me(), d, work, n, mpitype, target) \
               << endl; \
        } \
      ret = MPI_Reduce(d, work, n, mpitype, op, target, commgrp); \
      if (debug_) \
        ExEnv::outn() << scprintf("%3d: done with Reduce", me()) << endl; \
    } \
 \
  if (ret != MPI_SUCCESS) { \
      ExEnv::outn() << me() << ": MPIMessageGrp::reduce(," \
          << n << ",,," << target << "): mpi error:" << endl; \
      print_error_and_abort(me(), ret); \
    } \
 \
  if (target == -1 || target == me()) { \
     for (int i=0; i<n; i++) d[i] = work[i]; \
    } \
 \
  MPI_Op_free(&op); \
 \
  if (!scratch) delete[] work; \
}

REDUCEMEMBER(double, double, MPI_DOUBLE)
REDUCEMEMBER(float, float, MPI_FLOAT)
REDUCEMEMBER(uint, unsigned int, MPI_INT)
REDUCEMEMBER(int, int, MPI_INT)
REDUCEMEMBER(short, short, MPI_SHORT)
REDUCEMEMBER(long, long, MPI_LONG)
REDUCEMEMBER(char, char, MPI_CHAR)
REDUCEMEMBER(uchar, unsigned char, MPI_UNSIGNED_CHAR)
#ifdef MPI_SIGNED_CHAR
REDUCEMEMBER(schar, signed char, MPI_SIGNED_CHAR)
#else
void
MPIMessageGrp::reduce(signed char* d, int n, GrpReduce<signed char>& r,
                      signed char*scratch, int target)
{
  MessageGrp::reduce(d,n,r,scratch,target);
}
#endif

#define SUMMEMBER(name, type, mpitype) \
void \
MPIMessageGrp::sum(type*d, int n, type*scratch, int target) \
{ \
  if (use_messagegrp_collectives_) { \
      MessageGrp::sum(d,n,scratch,target); \
      return; \
    } \
 \
  type *work; \
  if (!scratch) work = new type[n]; \
  else work = scratch; \
 \
  int ret; \
 \
  if (target == -1) { \
      if (debug_) { \
          ExEnv::outn() << scprintf("%3d: MPI_Allreduce" \
          "(0x%08x, 0x%08x, %5d, %3d, MPI_SUM, commgrp)", \
          me(), d, work, n, mpitype) \
               << endl; \
        } \
      ret = MPI_Allreduce(d, work, n, mpitype, MPI_SUM, commgrp); \
      if (debug_) \
        ExEnv::outn() << scprintf("%3d: done with Allreduce", me()) << endl; \
    } \
  else { \
      if (debug_) { \
          ExEnv::outn() << scprintf("%3d: MPI_Reduce" \
          "(0x%08x, 0x%08x, %5d, %3d, MPI_SUM, %3d, commgrp)", \
          me(), d, work, n, mpitype, target) \
               << endl; \
        } \
      ret = MPI_Reduce(d, work, n, mpitype, MPI_SUM, target, commgrp); \
      if (debug_) \
        ExEnv::outn() << scprintf("%3d: done with Reduce", me()) << endl; \
    } \
 \
  if (ret != MPI_SUCCESS) { \
      ExEnv::outn() << me() << ": MPIMessageGrp::sum(," \
          << n << ",,," << target << "): mpi error:" << endl; \
      print_error_and_abort(me(), ret); \
    } \
 \
  if (target == -1 || target == me()) { \
     for (int i=0; i<n; i++) d[i] = work[i]; \
    } \
 \
  if (!scratch) delete[] work; \
}
SUMMEMBER(double, double, MPI_DOUBLE)
SUMMEMBER(int, int, MPI_INT)

void
MPIMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  if (n() == 1) return;

  if (use_messagegrp_collectives_) {
      MessageGrp::raw_bcast(data,nbyte,from);
      return;
    }

  if (debug_) {
      ExEnv::outn() << scprintf("%3d: MPI_Bcast("
                       "0x%08x, %5d, MPI_BYTE, %3d, commgrp)",
                       me(), data, nbyte, from)
           << endl;
    }
  int ret;
  if ((ret = MPI_Bcast(data, nbyte, MPI_BYTE, from, commgrp))
      != MPI_SUCCESS) {
      ExEnv::outn() << me() << ": MPIMessageGrp::raw_bcast(,"
          << nbyte << "," << from << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) {
      ExEnv::outn() << scprintf("%3d: done with bcast", me()) << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
