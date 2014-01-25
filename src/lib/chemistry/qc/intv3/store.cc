//
// store.cc
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

#include <stdlib.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/types.h>

#include <chemistry/qc/intv3/storage.h>
#include <chemistry/qc/intv3/int2e.h>

#define PRINT_STORED 0
#define MONITOR_HASH -1

using namespace sc;

void
Int2eV3::init_storage(int size)
{
  storer = new IntegralStorer();
  storer->init(size);
  if (size) int_integral_storage = size;
  else int_integral_storage = 0;
}

void
Int2eV3::done_storage()
{
  if (storer.nonnull()) {
      storer->done();
    }
  int_integral_storage = 0;
}

int
Int2eV3::int_have_stored_integral(int sh1,int sh2,int sh3,int sh4,
                                  int p12,int p34,int p13p24)
{
  IntegralKey key(sh1,sh2,sh3,sh4,p12,p34,p13p24);
  IntegralLink *integral = storer->find(key);

  if (!integral) return 0;

#if PRINT_STORED
  if (sh1 != integral->intlist.key.sh0()
      ||sh2 != integral->intlist.key.sh1()
      ||sh3 != integral->intlist.key.sh2()
      ||sh4 != integral->intlist.key.sh3()
      ||p12 != integral->intlist.key.p12()
      ||p34 != integral->intlist.key.p34()
      ||p13p24 != integral->intlist.key.p13p24()) {
      ExEnv::outn() << scprintf("!!!!! SHELL INFO INCONSISTENCY\n");
      abort();
    }
  ExEnv::outn() << scprintf("===== %2d %2d %2d %2d, %d %d %d size %5d cost %7d slot %5d\n",
         sh1, sh2, sh3, sh4,
         p12, p34, p13p24,
         integral->size, integral->costlist.key,
         integral->hash()%storer->table_size());
  if (integral->hash()%storer->table_size() == MONITOR_HASH) {
      storer->table_entry(MONITOR_HASH).detailed_print();
    }
#endif

  int i;
  double *buffer = integral->buffer();
  for (i=0; i<integral->size; i++) {
      int_buffer[i] = buffer[i];
    }
  
  return 1;
}

void
Int2eV3::int_store_integral(int sh1,int sh2,int sh3,int sh4,
                            int p12,int p34,int p13p24,
                            int size)
{
  // the cost of the integral is the time to evaluate it
  // times the number of times it is needed
  // divided by the amount of memory required to store it
  int cost;
  if (int_Qvec) cost = erep_4bound(sh1,sh2,sh3,sh4) + 30;
  else cost = 1;
  if (cost <= 0) return;
  cost *=  int_shell1->nprimitive()
         * int_shell2->nprimitive()
         * int_shell3->nprimitive()
         * int_shell4->nprimitive()
         * size
         * 1024; // the 1024 is arbitrary
  int actualsize = IntegralLink::size_to_actualsize(size);
  cost /= actualsize;

  if (storer->should_store(cost, actualsize)) {
      IntegralKey key(sh1,sh2,sh3,sh4,p12,p34,p13p24);
      storer->store(key,int_buffer,size,cost,actualsize);
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
