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

#include <stdio.h> // for sprintf
#include <unistd.h> // for fchdir etc.
#include <fcntl.h> // for open on AIX

#include <mpi.h>
extern int MPI_Initialized(int *); // missing in mpi.h

#include <util/keyval/keyval.h>
#include <util/group/messmpi.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>

using namespace std;

//#define MPI_SEND_ROUTINE MPI_Ssend // hangs
#define MPI_SEND_ROUTINE MPI_Send // hangs in old MPI implementations
//#define MPI_SEND_ROUTINE MPI_Bsend // works requires the attach and detach
#define MPI_SEND_ROUTINE_NAME "MPI_Send"

// OP_COMMUTES is zero to work around a bug in MPI/Pro 1.5b5 and earlier
#define OP_COMMUTES 1

///////////////////////////////////////////////////////////////////////

static
void
print_error_and_abort(int me, int mpierror)
{
  char msg[MPI_MAX_ERROR_STRING+1];
  int size;
  MPI_Error_string(mpierror, msg, &size);
  msg[size] = '\0';
  ExEnv::out() << me << ": " << msg << endl;
  ExEnv::out().flush();
  //MPI_Abort(MPI_COMM_WORLD, mpierror);
}

///////////////////////////////////////////////////////////////////////
// The MPIMessageGrp class

static ClassDesc MPIMessageGrp_cd(
  typeid(MPIMessageGrp),"MPIMessageGrp",1,"public MessageGrp",
  create<MPIMessageGrp>, create<MPIMessageGrp>, 0);

MPIMessageGrp::MPIMessageGrp()
{
  init();
}

MPIMessageGrp::MPIMessageGrp(const Ref<KeyVal>& keyval):
  MessageGrp(keyval)
{
  int argc = -1;
  char **argv = 0;
  if (keyval->exists("argv")) {
      argc = keyval->count("argv");
      argv = new char*[argc+1];
      argv[argc] = 0;
      for (int arg=0; arg<argc; arg++) {
          argv[arg] = keyval->pcharvalue("argv",arg);
        }
    }

  init(argc, argv);

  if (keyval->booleanvalue("errors_return")) {
      if (me()==0)
          ExEnv::out() << indent << "MPIMessageGrp: errors_return is true" << endl;
      MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    }

  SCFormIO::init_mp(me());
  if (debug_) {
      ExEnv::out() << indent << "MPIMessageGrp: KeyVal CTOR: done" << endl;
    }
}

void
MPIMessageGrp::init(int argc,char **argv)
{
  int me, nproc;

  if (debug_) {
      ExEnv::out() << "MPIMessageGrp::init: entered" << endl;
    }

  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
      if (argc < 0) {
          argc = 1;
          argv = new char*[argc+1];
          // reduce the internal buffer since a user buffer is used
          //argv[0] = "-mpiB4";
          //argv[1] = 0;
          argc = 0;
          argv[0] = 0;
        }
      // This dot business is to work around problems with some MPI
      // implementations.
      int dot = open(".",O_RDONLY);
      if (debug_) {
          ExEnv::out() << indent
               << "Calling MPI_Init with";
          for (int i=0; i<argc; i++) {
              ExEnv::out() << " " << argv[i];
            }
          ExEnv::out() << endl;
        }
      MPI_Init(&argc, &argv);
#ifdef HAVE_FCHDIR
      fchdir(dot);
#endif
      close(dot);
    }
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  bufsize = 4000000;
  buf = 0;
  //buf = (void*) new char[bufsize];
  //MPI_Buffer_attach(buf,bufsize);
  
  initialize(me, nproc);

  //MPIL_Trace_on();

  if (debug_) {
      ExEnv::out() << me << ": MPIMessageGrp::init: done" << endl;
    }
}

MPIMessageGrp::~MPIMessageGrp()
{
  //MPIL_Trace_off();
  //MPI_Buffer_detach(&buf, &bufsize);
  delete[] (char*) buf;
  MPI_Finalize();
}

void
MPIMessageGrp::raw_send(int target, void* data, int nbyte)
{
  if (debug_) {
      ExEnv::out() << scprintf("%3d: " MPI_SEND_ROUTINE_NAME
                       "(0x%08x, %5d, MPI_BYTE, %3d, 0, MPI_COMM_WORLD)",
                       me(), data, nbyte, target)
           << endl;
    }
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,0,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::raw_send("
          << target << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) ExEnv::out() << scprintf("%3d: sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  MPI_Status status;
  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: MPI_Recv"
                       "(0x%08x, %5d, MPI_BYTE, %3d, 0, MPI_COMM_WORLD,)",
                       me(), data, nbyte, sender)
           << endl;
    }
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,sender,0,MPI_COMM_WORLD,&status))
      != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::raw_recv("
          << sender << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) ExEnv::out() << scprintf("%3d: recvd %d bytes\n", me(), rlen) << endl;
}

void
MPIMessageGrp::raw_sendt(int target, int type, void* data, int nbyte)
{
  type = (type<<1) + 1;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: " MPI_SEND_ROUTINE_NAME
                       "(0x%08x, %5d, MPI_BYTE, %3d, %5d, MPI_COMM_WORLD)",
                       me(), data, nbyte, target, type)
           << endl;
    }
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,type,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::raw_sendt("
          << target << "," << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) ExEnv::out() << scprintf("%3d: sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  MPI_Status status;
  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: MPI_Recv(0x%08x, %5d, MPI_BYTE, "
                       "MPI_ANY_SOURCE, %5d, MPI_COMM_WORLD,)",
                       me(), data, nbyte, type)
           << endl;
    }
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,MPI_ANY_SOURCE,
                      type,MPI_COMM_WORLD,&status)) != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::raw_recvt("
          << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: recvd %d bytes from %d with tag %d\n",
                       me(), rlen, rnode, rtag) << endl;
    }
}

int
MPIMessageGrp::probet(int type)
{
  int flag;
  MPI_Status status;

  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  int ret;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: MPI_Iprobe(MPI_ANY_SOURCE, %5d, MPI_COMM_WORLD, "
                       "&flag, &status)", me(), type)
           << endl;
    }
  if ((ret = MPI_Iprobe(MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&flag,&status))
      != MPI_SUCCESS ) {
      ExEnv::out() << me() << ": MPIMessageGrp::probet("
          << type << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (flag) {
    rnode = status.MPI_SOURCE;
    rtag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_BYTE, &rlen);
    return 1;
    }
  else {
    rnode = rtag = rlen = 0;
    }
    
  return 0;
}

void
MPIMessageGrp::sync()
{
  int ret;
  if (debug_) {
      ExEnv::out() << scprintf("%3d: MPI_Barrier(MPI_COMM_WORLD)", me()) << endl;
    }
  if ((ret = MPI_Barrier(MPI_COMM_WORLD)) != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::sync(): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
}

#define REDUCEMEMBER(name, type, mpitype) \
static GrpReduce<type>* name ## reduceobject; \
static void \
name ## reduce(void*b, void*a, int*len, MPI_Datatype*datatype) \
{ \
  name ## reduceobject->reduce((type*)a, (type*)b, *len); \
} \
void \
MPIMessageGrp::reduce(type*d, int n, GrpReduce<type>&r, \
                      type*scratch, int target) \
{ \
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
          ExEnv::out() << scprintf("%3d: MPI_Allreduce" \
          "(0x%08x, 0x%08x, %5d, %3d, op, MPI_COMM_WORLD)", \
          me(), d, work, n, mpitype) \
               << endl; \
        } \
      ret = MPI_Allreduce(d, work, n, mpitype, op, MPI_COMM_WORLD); \
      if (debug_) \
        ExEnv::out() << scprintf("%3d: done with Allreduce", me()) << endl; \
    } \
  else { \
      if (debug_) { \
          ExEnv::out() << scprintf("%3d: MPI_Reduce" \
          "(0x%08x, 0x%08x, %5d, %3d, op, %3d, MPI_COMM_WORLD)", \
          me(), d, work, n, mpitype, target) \
               << endl; \
        } \
      ret = MPI_Reduce(d, work, n, mpitype, op, target, MPI_COMM_WORLD); \
      if (debug_) \
        ExEnv::out() << scprintf("%3d: done with Reduce", me()) << endl; \
    } \
 \
  if (ret != MPI_SUCCESS) { \
      ExEnv::out() << me() << ": MPIMessageGrp::reduce(," \
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

void
MPIMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  if (n() == 1) return;

  if (debug_) {
      ExEnv::out() << scprintf("%3d: MPI_Bcast("
                       "0x%08x, %5d, MPI_BYTE, %3d, MPI_COMM_WORLD)",
                       me(), data, nbyte, from)
           << endl;
    }
  int ret;
  if ((ret = MPI_Bcast(data, nbyte, MPI_BYTE, from, MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      ExEnv::out() << me() << ": MPIMessageGrp::raw_bcast(,"
          << nbyte << "," << from << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) {
      ExEnv::out() << scprintf("%3d: done with bcast", me()) << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
