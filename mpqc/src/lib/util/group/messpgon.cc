//
// messpgon.cc
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

extern "C" {
#include <nx.h>
}
#include <util/group/messpgon.h>

typedef long(*gopfarg)(...);

static ClassDesc ParagonMessageGrp_cd(
  typeid(ParagonMessageGrp),"ParagonMessageGrp",1,"public intMessageGrp",
  0, create<ParagonMessageGrp>, 0);

ParagonMessageGrp::ParagonMessageGrp(const Ref<KeyVal>&)
{
  // an extra control bit is needed for the bcast type
  ctl_nbit = 3;
  ctl_mask = 0x7;
  initialize();
}

ParagonMessageGrp::ParagonMessageGrp()
{
  // an extra control bit is needed for the bcast type
  ctl_nbit = 3;
  ctl_mask = 0x7;
  initialize();
}

void
ParagonMessageGrp::initialize()
{
  int nprocs = numnodes();
  int mynodeid = mynode();
  intMessageGrp::initialize(mynodeid, nprocs, 30);
}

int
ParagonMessageGrp::basic_probe(int type)
{
  return iprobe(type);
}

ParagonMessageGrp::~ParagonMessageGrp()
{
}

void
ParagonMessageGrp::sync()
{
  gsync();
}

void
ParagonMessageGrp::basic_recv(int type, void* buf, int bytes)
{
  crecv(type, (char*) buf, bytes);
}

void
ParagonMessageGrp::basic_send(int dest, int type, void* buf, int bytes)
{
  csend(type, (char*) buf, bytes, dest, 0);
}
 
int
ParagonMessageGrp::last_source()
{
  return infonode();
}

int
ParagonMessageGrp::last_size()
{
  return infocount();
}

int
ParagonMessageGrp::last_type()
{
  return msgtype_typ(infotype());
}

static int nreduce;

static GrpReduce<double>* doublereduceobject;
static long
doublereduce(char*a, char*b)
{
  doublereduceobject->reduce((double*)a, (double*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(double*d, int n, GrpReduce<double>&r,
                          double*scratch, int target)
{
  nreduce = n;
  doublereduceobject = &r;

  double *work;
  if (!scratch) work = new double[n];
  else work = scratch;

  gopf(d, n*sizeof(double), work, doublereduce);

  if (!scratch) delete[] work;
}

static GrpReduce<int>* intreduceobject;
static long
intreduce(char*a, char*b)
{
  intreduceobject->reduce((int*)a, (int*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(int*d, int n, GrpReduce<int>&r,
                          int*scratch, int target)
{
  nreduce = n;
  intreduceobject = &r;

  int *work;
  if (!scratch) work = new int[n];
  else work = scratch;

  gopf(d, n*sizeof(int), work, intreduce);

  if (!scratch) delete[] work;
}

static GrpReduce<char>* charreduceobject;
static long
charreduce(char*a, char*b)
{
  charreduceobject->reduce((char*)a, (char*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(char*d, int n, GrpReduce<char>&r,
                          char*scratch, int target)
{
  nreduce = n;
  charreduceobject = &r;

  char *work;
  if (!scratch) work = new char[n];
  else work = scratch;

  gopf(d, n*sizeof(char), work, charreduce);

  if (!scratch) delete[] work;
}

static GrpReduce<long>* longreduceobject;
static long
longreduce(char*a, char*b)
{
  longreduceobject->reduce((long*)a, (long*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(long*d, int n, GrpReduce<long>&r,
                          long*scratch, int target)
{
  nreduce = n;
  longreduceobject = &r;

  long *work;
  if (!scratch) work = new long[n];
  else work = scratch;

  gopf(d, n*sizeof(long), work, longreduce);

  if (!scratch) delete[] work;
}

static GrpReduce<float>* floatreduceobject;
static long
floatreduce(char*a, char*b)
{
  floatreduceobject->reduce((float*)a, (float*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(float*d, int n, GrpReduce<float>&r,
                          float*scratch, int target)
{
  nreduce = n;
  floatreduceobject = &r;

  float *work;
  if (!scratch) work = new float[n];
  else work = scratch;

  gopf(d, n*sizeof(float), work, floatreduce);

  if (!scratch) delete[] work;
}

static GrpReduce<unsigned char>* unsignedcharreduceobject;
static long
unsignedcharreduce(char*a, char*b)
{
  unsignedcharreduceobject->reduce((unsigned char*)a, (unsigned char*)b,
                                   nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(unsigned char*d, int n, GrpReduce<unsigned char>&r,
                          unsigned char*scratch, int target)
{
  nreduce = n;
  unsignedcharreduceobject = &r;

  unsigned char *work;
  if (!scratch) work = new unsigned char[n];
  else work = scratch;

  gopf(d, n*sizeof(unsigned char), work, unsignedcharreduce);

  if (!scratch) delete[] work;
}

static GrpReduce<short>* shortreduceobject;
static long
shortreduce(char*a, char*b)
{
  shortreduceobject->reduce((short*)a, (short*)b, nreduce);
  return 0;
}
void
ParagonMessageGrp::reduce(short*d, int n, GrpReduce<short>&r,
                          short*scratch, int target)
{
  nreduce = n;
  shortreduceobject = &r;

  short *work;
  if (!scratch) work = new short[n];
  else work = scratch;

  gopf(d, n*sizeof(short), work, shortreduce);

  if (!scratch) delete[] work;
}

void
ParagonMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  if (numnodes() == 1) return;

  static int bcast_type = 0;

  bcast_type++;
  if (bcast_type > typ_mask) bcast_type = 0;

  // find an unused message type
  int type = bcast_type<<typ_shift | 4<<ctl_shift;
  if (from == mynode()) {
      csend(type,(char*)data,nbyte,-1,0);
    }
  else {
      crecv(type,(char*)data,nbyte);
    }
}

void
ParagonMessageGrp::raw_collect(const void *part, const int *lengths, void *whole,
                              int bytes_per_datum)
{
  if (bytes_per_datum != 1) {
      int *newlengths = new int[n()];
      for (int i=0; i<n(); i++) newlengths[i] = lengths[i] * bytes_per_datum;
      gcolx((void*)part,(long*)newlengths,whole);
      delete[] newlengths;
    }
  else {
      gcolx((void*)part,(long*)lengths,whole);
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
