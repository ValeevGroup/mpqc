//
// memmpi.cc
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

#ifndef _util_group_memmpi_cc
#define _util_group_memmpi_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmpi.h>

#include <mpi.h>

///////////////////////////////////////////////////////////////////////
// The MPIMemoryGrp class

#define CLASSNAME MPIMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
MPIMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPIMemoryGrp::MPIMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  if (debug_) cout << "MPIMemoryGrp entered" << endl;

  use_acknowledgments_ = 0;
  use_active_messages_ = 0;

  init_mid();

  activate();
}

MPIMemoryGrp::MPIMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
{
  if (debug_) cout << "MPIMemoryGrp keyval entered" << endl;

  use_acknowledgments_ = 0;
  use_active_messages_ = 0;

  init_mid();

  activate();
}

void
MPIMemoryGrp::init_mid()
{
  for (int i=0; i<max_mid; i++) mid_ready_[i] = 1;
}

long
MPIMemoryGrp::get_mid()
{
  for (int i=0; i<max_mid; i++) {
      if (mid_ready_[i]) {
          mid_ready_[i] = 0;
          if (debug_) cout << "MPIMemoryGrp::get_mid(): got " << i << endl;
          return i;
        }
    }

  cerr << "MPIMemoryGrp::get_mid(): ran out of mid's" << endl;
  abort();
  return 0;
}

void
MPIMemoryGrp::free_mid(long mid)
{
  if (debug_) cout << "MPIMemoryGrp::free_mid(): freeing " << mid << endl;
  mid_ready_[mid] = 1;
}

long
MPIMemoryGrp::lockcomm()
{
  return 0;
}

void
MPIMemoryGrp::unlockcomm(long oldvalue)
{
}

long
MPIMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  int mid = get_mid();
  //int MPI_Ibsend(void* buf, int count, MPI_Datatype datatype,
  //int dest, int tag, MPI_Comm comm, MPI_Request *request) 
  if (debug_) cout << ">>>> MPI_Ibsend for mid " << mid
                   << " type " << type << endl;
  MPI_Ibsend(data, nbytes, MPI_BYTE, node, type,
             MPI_COMM_WORLD, &handles_[mid]);
  return mid;
}

long
MPIMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  int n;
  if (node == -1) n = MPI_ANY_SOURCE;
  else n = node;
  int t = type;
  int mid = get_mid();
  // int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int
  // source, int tag, MPI_Comm comm, MPI_Request *request) 
  if (debug_) cout << ">>>> MPI_Irecv for mid " << mid
                   << " type " << type << endl;
  MPI_Irecv(data, nbytes, MPI_BYTE, n, t,
            MPI_COMM_WORLD, &handles_[mid]);
  if (debug_) cerr << "MPIMemoryGrp:: recv mid = " << mid << endl;
  return mid;
}

long
MPIMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  cerr << "MPIMemoryGrp::postrecv: active messages not supported\n" << endl;
  abort();
  return 0;
}

long
MPIMemoryGrp::wait(long mid1, long mid2)
{
  MPI_Status status;
  if (debug_) {
      cout << me() << ": MPIMemoryGrp::wait("
           << mid1 << ", "
           << mid2 << ")"
           << endl;
    }

  if (mid1 == -1) {
      cerr << me() << ": MPIMemoryGrp::wait: mid1 == -1" << endl;
      sleep(1);
      abort();
    }
  else if (mid2 == -1) {
      if (debug_) cout << ">>>> MPI_Wait for " << mid1 << endl;
      MPI_Wait(&handles_[mid1], &status);
      free_mid(mid1);
      if (debug_)
          cout << me() << ": MPIMemoryGrp::wait(): got(1) " << mid1 << endl;
      return mid1;
    }
  else {
      while(1) {
          int flag;
          //if (debug_) cout << ">>>> MPI_Test for " << mid1 << endl;
          MPI_Test(&handles_[mid1], &flag, &status);
          if (flag) {
              free_mid(mid1);
              if (debug_)
                  cout << me() << ": MPIMemoryGrp::wait(): got(2a) "
                       << mid1 << endl;
              return mid1;
            }
          //if (debug_) cout << ">>>> MPI_Test for " << mid2 << endl;
          MPI_Test(&handles_[mid2], &flag, &status);
          if (flag) {
              free_mid(mid2);
              if (debug_)
                  cout << me() << ": MPIMemoryGrp::wait(): got(2b) "
                       << mid2 << endl;
              return mid2;
            }
        }
    }
}

int
MPIMemoryGrp::probe(long mid)
{
  MPI_Status status;
  int flag;
  MPI_Test(&handles_[mid], &flag, &status);
  if (flag) {
      free_mid(mid);
      if (debug_)
          cout << me() << ": MPIMemoryGrp::probe(): got "
               << mid << endl;
      return 1;
    }
  return 0;
}

void
MPIMemoryGrp::deactivate()
{
  if (active_) {
      if (debug_) cout << "MPIMemoryGrp::deactivate()" << endl;
      sync();
      if (debug_) cout << ">>>> MPI_Cancel for " << data_request_mid_ << endl;
      MPI_Cancel(&handles_[data_request_mid_]);
      free_mid(data_request_mid_);
      active_ = 0;
    }
}

MPIMemoryGrp::~MPIMemoryGrp()
{
  if (debug_) cerr << "MPIMemoryGrp: in DTOR" << endl;
  deactivate();
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
