
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_types.h>

#include <chemistry/qc/intv2/storage.h>
#include <chemistry/qc/intv2/inter.h>
#include <chemistry/qc/intv2/bounds.h>

#include <chemistry/qc/intv2/bounds.gbl>

#define PRINT_STORED 0

/////////////////////////////////////////////////////////////////////////
// C language interface to a global integral storer

static RefIntegralStorer storer;

void
int_storage(int size)
{
  storer = new IntegralStorer();
  storer->init(size);
  if (size) int_integral_storage = 1;
  else int_integral_storage = 0;
}

void
int_reduce_storage_threshold()
{
}

void
int_done_storage()
{
  if (storer.nonnull()) {
      storer->done();
    }
  int_integral_storage = 0;
}

int
int_have_stored_integral(int sh1,int sh2,int sh3,int sh4,
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
      printf("!!!!! SHELL INFO INCONSISTENCY\n");
      abort();
    }
  printf("===== %d %d %d %d, %d %d %d size %5d cost %7d at 0x%x slot %5d\n",
         sh1, sh2, sh3, sh4,
         p12, p34, p13p24,
         integral->size, integral->costlist.key,
         integral, integral->hash()%storer->table_size());
#endif

  int i;
  double *buffer = integral->buffer();
  for (i=0; i<integral->size; i++) {
      int_buffer[i] = buffer[i];
    }
  
  return 1;
}

void
int_store_integral(int sh1,int sh2,int sh3,int sh4,
                   int p12,int p34,int p13p24,
                   int size)
{
  // the cost of the integral is the time to evaluate it
  // times the number of times it is needed
  // divided by the amount of memory required to store it
  int cost;
  if (int_Qvec) cost = int_erep_4bound(sh1,sh2,sh3,sh4) + 30;
  else cost = 1;
  if (cost <= 0) return;
  cost *=  int_shell1->nprim
         * int_shell2->nprim
         * int_shell3->nprim
         * int_shell4->nprim
         * size
         * 1024; // the 1024 is arbitrary
  int actualsize = IntegralLink::size_to_actualsize(size);
  cost /= actualsize;

  if (storer->should_store(cost, actualsize)) {
      IntegralKey key(sh1,sh2,sh3,sh4,p12,p34,p13p24);
      storer->store(key,int_buffer,size,cost,actualsize);
    }
}
