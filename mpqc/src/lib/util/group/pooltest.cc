
#include <stdlib.h>
#include <math.h>
#include <ostream.h>
#include <util/group/pool.h>

class Double {
  private:
    static Pool* pool_;
    static Double* list;
    Double* next;
    Double* prev;
    double* d;
    int size;
    static void zaplist(Double*);
  public:
    Double(size_t size);
    ~Double();
    void zap();
    static void zapall();
    void clear();
    static void pool(Pool*);
};

Double* Double::list = 0;
Pool* Double::pool_ = 0;

void
Double::clear()
{
  if (d) pool_->release(d);
  d = 0;
  size = 0;
}

void
Double::zapall()
{
  zaplist(list);
}

void
Double::zaplist(Double*l)
{
  for (Double* i=l; i; i = i->next) {
      i->zap();
    }
}

Double::Double(size_t s):
  size(s)
{
  if (!pool_) {
      fprintf(stderr,"Double::Double: Pool not initialized\n");
      abort();
    }
  d = pool_->allocate_double(size);
  if (!d) {
      //fprintf(stdout,"\nDouble::Double allocation of size %d failed\n",size);
      fprintf(stdout,"F", size); fflush(stdout);
      size = 0;
    }
  next = list;
  prev = 0;
  list = this;
  if (next) next->prev = this;
}

Double::~Double()
{
  clear();
  if (next) next->prev = prev;
  if (prev) prev->next = next;
  else list = next;
}

void
Double::zap()
{
  int* x = (int*)d;
  for (int i=0; i<size*2; i++) {
      if (x[i] == PoolData::magic) {
          fprintf(stderr,"Double::zap: tried to zap a magic number\n");
          abort();
        }
    }
  for (i=0; i<size; i++) d[i] = 0.0;
}

void
Double::pool(Pool*p)
{
  if (pool_ && list) {
      fprintf(stderr,"Double::pool: cannot reinitialize pool\n");
      abort();
    }
  pool_ = p;
}

void test1(Pool*);
void test2(Pool*);

int
main()
{
  const int poolsize = 4000000;
  Pool *pool = new(malloc(poolsize)) Pool(poolsize);

  Double::pool(pool);

  srand48(100);

  printf("test1:\n");
  test1(pool);
  printf("test2:\n");
  test2(pool);

  return 0;
}


void
test1(Pool*pool)
{
  pool->check();

  pool->print();
  printf("\n");

  Double t1(10);

  Double::zapall();

  pool->check();

  pool->print();
  printf("\n");

  Double t2(10000);

  Double::zapall();

  pool->check();

  Double t3(100);

  Double::zapall();

  pool->check();
  
  pool->print();
  printf("\n");

  Double::zapall();

  pool->check();
  
  pool->print();
  printf("\n");

  t2.clear();

  pool->check();
  
  pool->print();
  printf("\n");

  Double t4(100);

  pool->check();

  pool->print();
  printf("\n");

  t1.clear();
  t4.clear();
  t3.clear();

  pool->check();

  pool->print();
  printf("\n");

}

void
test2(Pool*pool)
{
  const int npass = 200;
  const int nd = 4096;
  Double* d[nd];

  for (int i=0; i<nd; i++) d[i] = 0;

  for (int ii=0; ii<npass; ii++) {
      for (i=0; i<nd;) {
          if (mrand48() > 0) {
              // allocate data
              size_t size = lrand48() & 0x03ff;
              d[i] = new Double(size);
              fflush(stderr); printf("a"); fflush(stdout);
              i++;
            }
          else {
              // deallocate data
              int loc = (int) (drand48()*i);
              if (loc >= nd) loc = nd-1;
              if (loc < 0) loc = 0;
              if (d[loc]) {
                  fflush(stderr); printf("d"); fflush(stdout);
                  delete d[loc];
                  d[loc] = 0;
                }
            }
          //pool->print();
          //pool->check();
          //Double::zapall();
        }
      for (i=0; i<nd; i++) {
          if (d[i]) {
              fflush(stderr); printf("d"); fflush(stdout);
              delete d[i];
              d[i] = 0;
            }
        }
      printf("\n");
      pool->print();
      //pool->check();
      //Double::zapall();
    }
}
