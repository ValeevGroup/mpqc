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

#include <mpi.h>
#include <util/keyval/keyval.h>
#include <util/group/messmpi.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>

#if HAVE_P4
#include <netdb.h>
#endif

//#define MPI_SEND_ROUTINE MPI_Ssend // hangs
//#define MPI_SEND_ROUTINE MPI_Send // hangs
#define MPI_SEND_ROUTINE MPI_Bsend // works requires the attach and detach

///////////////////////////////////////////////////////////////////////

static
void
print_error_and_abort(int me, int mpierror)
{
  char msg[MPI_MAX_ERROR_STRING+1];
  int size;
  MPI_Error_string(mpierror, msg, &size);
  msg[size] = '\0';
  cerr << me << ": " << msg << endl;
  cout.flush();
  cerr.flush();
  MPI_Abort(MPI_COMM_WORLD, mpierror);
  abort();
}

///////////////////////////////////////////////////////////////////////
// The MPIMessageGrp class

#define CLASSNAME MPIMessageGrp
#define PARENTS public MessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
MPIMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPIMessageGrp::MPIMessageGrp()
{
  init();
}

MPIMessageGrp::MPIMessageGrp(const RefKeyVal& keyval):
  MessageGrp(keyval)
{
#if HAVE_P4
  nlocal=keyval->intvalue("nlocal");
  if (!nlocal || keyval->error() != KeyVal::OK)
      nlocal=1;

  nremote=keyval->count("remote_clusters");
  if (nremote) {
      remote_clusters = new p4_cluster[nremote];
      for (int i=0; i < nremote; i++) {
          remote_clusters[i].hostname =
                                   keyval->pcharvalue("remote_clusters",i,0);
          if (keyval->error() != KeyVal::OK) {
              cerr << indent
                   << "MPIMessageGrp(RefKeyVal): bad hostname" << endl;
              abort();
            }

          remote_clusters[i].nslaves = keyval->intvalue("remote_clusters",i,1);
          if (keyval->error() != KeyVal::OK) {
              cerr << indent
                   << "MPIMessageGrp(RefKeyVal): bad nslaves" << endl;
              abort();
            }
        }
    }
  else {
      remote_clusters=0;
    }

  master=keyval->pcharvalue("master");
  if (keyval->error() != 0) {
      if (nremote) {
          cerr << indent
               << "MPIMessageGrp(RefKeyVal): need master" << endl;
          abort();
        }
      else {
          char hostname[80];
          if (gethostname(hostname,79) < 0) {
              cerr << indent
                   << "MPIMessageGrp(RefKeyVal): could not get hostname"
                   << endl;
              abort();
            }
          master=new_string(hostname);
        }
    }

  jobid=keyval->pcharvalue("job_id");
  if (keyval->error() != 0) {
      cerr << "MPIMessageGrp(RefKeyVal): no job_id in input" << endl;
      abort();
    }
#endif

  init();

  if (keyval->booleanvalue("errors_return")) {
      if (me()==0)
          cout << indent << "MPIMessageGrp: errors_return is true" << endl;
      MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    }

  SCFormIO::init_mp(me());
  if (debug_) {
      cerr << indent << "MPIMessageGrp: KeyVal CTOR: done" << endl;
    }
}

#if HAVE_P4

static struct hostent *
mygethostbyname(const char *mach)
{
  int i,j;
  char *alias;
  struct hostent *host;

  if (mach) {
      host=gethostbyname(mach);
    }
  else {
      char name[80];
      gethostname(name,80);
      host=gethostbyname(name);
    }

  if (!host)
      return 0;
  
  struct hostent *mhost = new hostent;

  mhost->h_addrtype=host->h_addrtype;
  mhost->h_length=host->h_length;
  mhost->h_name=new_string(host->h_name);

  if (host->h_aliases) {
      i=0;
      while (alias=host->h_aliases[i]) i++;
      mhost->h_aliases = new char*[i+1];
      for (j=0; j < i ; j++)
          mhost->h_aliases[j]=new_string(host->h_aliases[j]);
      mhost->h_aliases[i]=NULL;
    }

  return mhost;
}

static void
freehostent(struct hostent *h)
{
  if (!h) return;
  
  delete[] h->h_name;

  int i=0;
  while (h->h_aliases[i]) {
      delete[] h->h_aliases[i];
      i++;
    }

  delete[] h->h_aliases;
  delete h;
}

static int
same_host(const char *m1, const char *m2)
{
  int i;
  char *alias;
  struct hostent *h1, *h2;

  if (!strcmp(m1,"any") || !strcmp(m2,"any"))
      return 1;

  h1 = mygethostbyname(m1);
  h2 = mygethostbyname(m2);

  if (!h1 || !h2) {
      freehostent(h1);
      freehostent(h2);
      return 0;
    }
  
  if (!strcmp(h1->h_name,h2->h_name)) {
      freehostent(h1);
      freehostent(h2);
      return 1;
    }

  i=0;
  while (alias=h1->h_aliases[i]) {
      if (!strcmp(alias,h2->h_name)) {
          freehostent(h1);
          freehostent(h2);
          return 1;
        }
      i++;
    }

  i=0;
  while (alias=h2->h_aliases[i]) {
      if (!strcmp(alias,h1->h_name)) {
          freehostent(h1);
          freehostent(h2);
          return 1;
        }
      i++;
    }

  freehostent(h1);
  freehostent(h2);
  return 0;
}

struct MPIMessageGrp::p4_cluster *
MPIMessageGrp::my_node_info(const char hostname[], int& nodenum)
{
  for (nodenum=1; nodenum <= nremote; nodenum++) {
      if (same_host(hostname,remote_clusters[nodenum-1].hostname))
          return &remote_clusters[nodenum-1];
    }
  return 0;
}
#endif

void
MPIMessageGrp::init()
{
  int me, nproc;
  int argc;
  char **argv;

#if HAVE_P4
  int am_slave = (getenv("SCSLAVE")!= 0);
  char hostname[80];
  if (gethostname(hostname,79) < 0) {
      cerr << indent
           << "MPIMessageGrp::init: could not get hostname" << endl;
      abort();
    }

  char numprocs[12], numnodes[12], nodenum[12];
  char **nslaves = new char*[nremote];
  
  if (!am_slave) {
      sprintf(numprocs,"%d",nlocal);
      sprintf(numnodes,"%d",nremote+1);
      sprintf(nodenum,"%d",0);

      argc = 16+4*nremote;
      argv = new char*[argc+1];
      argv[0] = "scprog";
      argv[1] = "-execer_id";
      argv[2] = "SC";
      argv[3] = "-master_host";
      argv[4] = master;
      argv[5] = "-my_hostname";
      argv[6] = hostname; 
      argv[7] = "-my_nodenum";
      argv[8] = nodenum;
      argv[9] = "-my_numprocs";
      argv[10] = numprocs;
      argv[11] = "-total_numnodes";
      argv[12] = numnodes;
      argv[13] = "-jobid";
      argv[14] = jobid;
      argv[15] = "-remote_info";

      int a=16;
      for (int i=0; i < nremote; i++) {
          nslaves[i] = new char[12];
          sprintf(nslaves[i],"%d",remote_clusters[i].nslaves);

          argv[a++] = remote_clusters[i].hostname;
          argv[a++] = "0";
          argv[a++] = nslaves[i];
          argv[a++] = "0";
        }
      
      argv[argc] = 0;
    }
  else {
      int mynum;
      struct p4_cluster *my_node = my_node_info(hostname, mynum);
      if (!my_node) {
          cerr << indent
               << "MPIMessageGrp:init: my_node is null" << endl;
          abort();
        }
              
      sprintf(numprocs,"%d",my_node->nslaves);
      sprintf(numnodes,"%d",nremote+1);
      sprintf(nodenum,"%d",mynum);
      
      argc = 15;
      argv = new char*[argc+1];
      argv[0] = "scprog";
      argv[1] = "-execer_id";
      argv[2] = "SC";
      argv[3] = "-master_host";
      argv[4] = master;
      argv[5] = "-my_hostname";
      argv[6] = hostname; 
      argv[7] = "-my_nodenum";
      argv[8] = nodenum;
      argv[9] = "-my_numprocs";
      argv[10] = numprocs;
      argv[11] = "-total_numnodes";
      argv[12] = numnodes;
      argv[13] = "-jobid";
      argv[14] = jobid;
      argv[15] = 0;
    }
#else  
  argc = 1;
  argv = new char*[argc+1];
  argv[0] = "-mpiB4"; // reduce the internal buffer since a user buffer is used
  argv[1] = 0;
#endif

  if (debug_) {
      cerr << "MPIMessageGrp::init: entered" << endl;
    }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  bufsize = 4000000;
  buf = (void*) new char[bufsize];
  MPI_Buffer_attach(buf,bufsize);
  initialize(me, nproc);

  if (debug_) {
      cerr << me << ": MPIMessageGrp::init: done" << endl;
    }

#if HAVE_P4
  if (!am_slave)
      for (int i=0; i < nremote; i++)
          delete[] nslaves[i];
  delete[] nslaves;
#endif  
}

MPIMessageGrp::~MPIMessageGrp()
{
  MPI_Buffer_detach(&buf, &bufsize);
  delete[] (char*) buf;
  MPI_Finalize();

#if HAVE_P4
  if (master) {
      delete[] master;
      master=0;
    }

  if (jobid) {
      delete[] jobid;
      jobid=0;
    }

  if (nremote && remote_clusters) {
      for (int i=0; i < nremote; i++)
          delete[] remote_clusters[i].hostname;
      delete[] remote_clusters;
      remote_clusters=0;
    }
#endif
}

void
MPIMessageGrp::raw_send(int target, void* data, int nbyte)
{
  if (debug_) {
      cerr << scprintf("Node %d sending %d bytes to %d with tag %d\n",
                       me(), nbyte, target, 0)
           << endl;
    }
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,0,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_send("
          << target << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) cerr << scprintf("Node %d sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  MPI_Status status;
  if (sender == -1) sender = MPI_ANY_SOURCE;
  if (debug_) {
      cerr << scprintf("Node %d recving %d bytes from %d with tag %d\n",
                       me(), nbyte, sender, 0)
           << endl;
    }
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,sender,0,MPI_COMM_WORLD,&status))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_recv("
          << sender << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) cerr << scprintf("Node %d recvd %d bytes\n", me(), rlen) << endl;
}

void
MPIMessageGrp::raw_sendt(int target, int type, void* data, int nbyte)
{
  type = (type<<1) + 1;
  if (debug_) {
      cerr << scprintf("Node %d sending %d bytes to %d with tag %d\n",
                       me(), nbyte, target, type)
           << endl;
    }
  int ret;
  if ((ret = MPI_SEND_ROUTINE(data,nbyte,MPI_BYTE,target,type,MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_sendt("
          << target << "," << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) cerr << scprintf("Node %d sent\n", me()) << endl;
}

void
MPIMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  MPI_Status status;
  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  if (debug_ ) {
      cerr << scprintf("Node %d recving %d bytes from %d with tag %d\n",
                       me(), nbyte, MPI_ANY_SOURCE, type)
           << endl;
    }
  int ret;
  if ((ret = MPI_Recv(data,nbyte,MPI_BYTE,MPI_ANY_SOURCE,
                      type,MPI_COMM_WORLD,&status)) != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_recvt("
          << type << ",," << nbyte << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  rnode = status.MPI_SOURCE;
  rtag = status.MPI_TAG;
  rlen = nbyte;
  if (debug_) cerr << scprintf("Node %d recvd %d bytes\n", me(), rlen) << endl;
}

int
MPIMessageGrp::probet(int type)
{
  int flag;
  MPI_Status status;

  if (type == -1) type = MPI_ANY_TAG;
  else type = (type<<1) + 1;
  int ret;
  if ((ret = MPI_Iprobe(MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&flag,&status))
      != MPI_SUCCESS ) {
      cerr << me() << ": MPIMessageGrp::probet("
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

int
MPIMessageGrp::last_source()
{
  return rnode;
}

int
MPIMessageGrp::last_size()
{
  return rlen;
}

int
MPIMessageGrp::last_type()
{
  return rtag>>1;
}

void
MPIMessageGrp::sync()
{
  int ret;
  if ((ret = MPI_Barrier(MPI_COMM_WORLD)) != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::sync(): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
}

static GrpReduce<double>* doublereduceobject;
static void
doublereduce(void*b, void*a, int*len, MPI_Datatype*datatype)
{
  doublereduceobject->reduce((double*)a, (double*)b, *len);
}
void
MPIMessageGrp::reduce(double*d, int n, GrpReduce<double>&r,
                          double*scratch, int target)
{
  doublereduceobject = &r;

  MPI_Op op;
  MPI_Op_create(doublereduce, 1, &op);

  double *work;
  if (!scratch) work = new double[n];
  else work = scratch;

  int ret;

  if (target == -1) {
      ret = MPI_Allreduce(d, work, n, MPI_DOUBLE, op, MPI_COMM_WORLD);
    }
  else {
      ret = MPI_Reduce(d, work, n, MPI_DOUBLE, op, target, MPI_COMM_WORLD);
    }

  if (ret != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::reduce(,"
          << n << ",,," << target << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }

  if (target == -1 || target == me()) {
     for (int i=0; i<n; i++) d[i] = work[i];
    }

  MPI_Op_free(&op);

  if (!scratch) delete[] work;
}


void
MPIMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  if (n() == 1) return;

  if (debug_) {
      cerr << scprintf("Node %d bcast %d bytes from %d", me(), nbyte, from)
           << endl;
    }
  int ret;
  if ((ret = MPI_Bcast(data, nbyte, MPI_BYTE, from, MPI_COMM_WORLD))
      != MPI_SUCCESS) {
      cerr << me() << ": MPIMessageGrp::raw_bcast(,"
          << nbyte << "," << from << "): mpi error:" << endl;
      print_error_and_abort(me(), ret);
    }
  if (debug_) {
      cerr << scprintf("Node %d done with bcast", me()) << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
